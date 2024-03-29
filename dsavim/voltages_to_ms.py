"""Utilties for converting T3 voltage dumps to measurement sets."""

from collections import namedtuple
from typing import Callable, Tuple
import time
import re
import os
from functools import wraps
from contextlib import contextmanager

import subprocess
from pkg_resources import resource_filename
import numpy as np
import astropy.units as u
import astropy.constants as c
from dask.distributed import Queue, Variable, Lock
from asyncio import TimeoutError
import asyncio

from antpos.utils import get_itrf
from dsautils.coordinates import get_declination, get_elevation
import dsacalib.constants as ct

from dsavim.dedisperse import get_dispersion_delay_ms
from dsavim.generate_uvh5 import calculate_uvw_and_geodelay, get_total_delay
from dsavim.utils import rsync_file, get_tstart_from_json, get_DM_from_json, load_params
from dsavim.generate_uvh5 import generate_uvh5

from dsamfs.fringestopping import generate_fringestopping_table

__all__ = [
    'pipeline_component', 'rsync_component', 'correlate_component',
    'uvh5_component',
    'generate_delay_table', 'initialize_system', 'initialize_candidate', 'initialize_correlator',
    'initialize_uvh5', 'initialize_vis_params', 'parse_visibility_parameters', 'get_cable_delays',
    'get_declination_etcd', 'get_blen', 'ProtectedVariable']

PARAMFILE = resource_filename('dsavim', 'data/voltage_corr_parameters.yaml')


def pipeline_component(targetfn, inqueue, outqueue=None):
    """Generate a component of the pipeline."""
    @wraps(targetfn)
    def inner():
        """Process data from a queue."""
        targetname = targetfn.__name__ if hasattr(targetfn, '__name__') else targetfn.func.__name__
        done = False
        while not done:
            # Get the next item
            try:
                item = inqueue.get(timeout=2)
                assert item
            except TimeoutError:
                time.sleep(8)
                continue

            # Process the item
            if item == 'END':
                done = True
            else:
                item = targetfn(item)

            # Pass the item on to the next queue
            if outqueue is not None:
                try:
                    print(f"{targetname} sending {item} to queue")
                    outqueue.put(item)
                except (EOFError, BrokenPipeError) as exc:
                    # Log and end gracefully if the queue is broken
                    print(
                        f"{type(exc).__name__} error when accessing outqueue "
                        f"in {targetname}")
                    done = True

    return inner


def rsync_component(item: Tuple[str], local: bool) -> str:
    """Rsync or copy a file."""
    srcfile, vfile = item
    if not os.path.exists(vfile):
        if not local:
            rsync_file(srcfile, vfile)
        else:
            os.symlink(srcfile, vfile)
    return vfile


def correlate_component(
        vfile: str, prev: str, dispersion_measure: float, ntint: int, corr_ch0: dict, npol: int, wfile: dict) -> str:

    """Correlate a file."""    
    sb = re.findall(r"sb\d\d", vfile)[0]
    _corrlist = {'sb00':'corr03','sb01':'corr04','sb02':'corr05','sb03':'corr06','sb04':'corr07',
                 'sb05':'corr08','sb06':'corr10','sb07':'corr11','sb08':'corr12','sb09':'corr14',
                 'sb10':'corr15','sb11':'corr16','sb12':'corr18','sb13':'corr19','sb14':'corr21',
                 'sb15':'corr22'}
    corr = _corrlist[sb]
    
    if not os.path.exists(f"{vfile}.corr"):
        first_channel_MHz = corr_ch0[corr]
        command = (
            "/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/toolkit "
            f"-i {vfile} -o {vfile}.corr -t {ntint} -c {first_channel_MHz} "
            f"-d delays.dat {'' if npol==4 else '-a'}")
        if dispersion_measure is not None:
            command += f" -m {dispersion_measure}"
        if wfile is not None:
            mywfile = f"{wfile['path']}/beamformer_weights_{sb}_{wfile['datestring']}.dat"
            command += f" -w {mywfile}"
        print(command)
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True)
        proc_stdout = str(process.communicate()[0].strip())
        print(proc_stdout)

    corrfile = f"{vfile}.corr"

    return corrfile


def uvh5_component(
        corrfile: str, candname: str, corrdir: str, declination: float,
        vis_params: dict, start_offset: int, end_offset: int) -> None:
    """Write correlated data to a uvh5 file."""
    _uvh5name = generate_uvh5(
        f"{corrdir}/{candname}",
        declination*u.deg,
        corrfile=corrfile,
        vis_params=vis_params,
        start_offset=start_offset,
        end_offset=end_offset
    )

    os.remove(corrfile)


def get_declination_etcd(tstart):
    """Look up the declination from etcd."""
    return get_declination(get_elevation(tstart)).to_value(u.deg)


def generate_delay_table(vis_params, reftime, declination, use_fs=True):
    """Generate a table of geometric and cable delays for the correlator."""

    if use_fs:
        # try using intermediate fringestopping table
        print("Using fringestopping table...")
        odelays = {"110": -222, "113": -1394, "114": -3514, "115": -5070}
        print(odelays)
        generate_fringestopping_table(vis_params['blen'],declination*np.pi/180.,1,1.,vis_params['antenna_order'],odelays,vis_params['bname'],reftime.mjd)
        data = np.load('fringestopping_table.npz')
        bws = -data['bw']/ct.C_GHZ_M 
        total_delay_string = '\n'.join(bws.flatten().astype('str'))+'\n' 

    else:
        _buvw, ant_bw = calculate_uvw_and_geodelay(vis_params, reftime, declination*u.deg)
        total_delay = get_total_delay(
            vis_params['baseline_cable_delays'], ant_bw, vis_params['bname'],
            vis_params['antenna_order'])
        total_delay_string = '\n'.join(total_delay.flatten().astype('str'))+'\n'
        
    with open('delays.dat', 'w', encoding='utf-8') as f:
        f.write(total_delay_string)


def initialize_system():
    """Set system parameters needed in voltage to ms converter."""
    params = load_params(PARAMFILE)
    SystemSetup = namedtuple(
        'SystemSetup', 'T3dir archivedir corrdir msdir start_time_offset reffreq_GHz corr_ch0_MHz')
    start_time_offset = params['burst_start_s']*u.s
    corr_ch0_MHz = {corr: 1e3*params['f0_GHz']+params['deltaf_MHz']*corr_ch0
                    for corr, corr_ch0 in params['ch0'].items()}

    system_setup = SystemSetup(
        params['T3dir'], params['archivedir'], params['corrdir'], params['msdir'],
        start_time_offset, params['reffreq_GHz'], corr_ch0_MHz)
    return system_setup


def initialize_candidate(candname, system_setup, dispersion_measure=None):
    """Set candidate parameters using information in the json header."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    sblist = ['sb00','sb01','sb02','sb03','sb04','sb05','sb06','sb07',
              'sb08','sb09','sb10','sb11','sb12','sb13','sb14','sb15']

    local = True
    headerfile = f"{system_setup.archivedir}/{candname}/Level2/voltages/T2_{candname}.json"
    voltagefiles = [f"{system_setup.archivedir}/{candname}/Level2/voltages/{candname}_{sb}_data.out"
                    for sb in sblist]

    tstart = get_tstart_from_json(headerfile) - system_setup.start_time_offset

    if dispersion_measure is None:
        dispersion_measure = get_DM_from_json(headerfile)
        if dispersion_measure < 1:
            dispersion_measure = None

    Candidate = namedtuple('Candidate', 'name time dm voltagefiles headerfile local')
    cand = Candidate(candname, tstart, dispersion_measure, voltagefiles, headerfile, local)
    return cand


def initialize_correlator(fullpol, ntint, cand, system_setup):
    """Set correlator parameters."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    sblist = ['sb00','sb01','sb02','sb03','sb04','sb05','sb06','sb07',
              'sb08','sb09','sb10','sb11','sb12','sb13','sb14','sb15']
    reftime = cand.time + system_setup.start_time_offset
    npol = 4 if fullpol else 2
    nfint = 1 if fullpol else 8
    corrfiles = [f"{system_setup.corrdir}/{cand.name}_{sb}_data.out" for sb in sblist]

    CorrelatorParameters = namedtuple('Correlator', 'reftime npol nfint ntint files')
    correlator_params = CorrelatorParameters(reftime, npol, nfint, ntint, corrfiles)
    return correlator_params


def initialize_uvh5(corrparams, cand, system_setup, outrigger_delays=None):
    """Set parameters for writing uvh5 files."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    sblist = ['sb00','sb01','sb02','sb03','sb04','sb05','sb06','sb07',
              'sb08','sb09','sb10','sb11','sb12','sb13','sb14','sb15']
    uvh5files = [f"{system_setup.corrdir}/{cand.name}_{sb}.hdf5" for sb in sblist]
    visparams = initialize_vis_params(corrparams, cand, outrigger_delays)

    UVH5Parameters = namedtuple('UVH5', 'files visparams')
    uvh5_params = UVH5Parameters(uvh5files, visparams)
    return uvh5_params


def initialize_vis_params(corrparams, cand, outrigger_delays=None):
    """Set parameters that describe the visbilities produced by the correlator.

    The time is adjusted so that the dedispersed signal has the correct time
    for the centre of the band.
    """
    T3params = load_params(PARAMFILE)
    if outrigger_delays:
        T3params['outrigger_delays'] = outrigger_delays

    # Adjust the start time so that the dedispersed signal has the correct
    # time for the centre of the band.
    dispersion_delay = get_dispersion_delay_ms(1.405, cand.dm, 1.530)*u.ms
    tstart = cand.time + dispersion_delay

    vis_params = parse_visibility_parameters(T3params, cand.time, corrparams.ntint)
    vis_params['tref'] = corrparams.reftime
    vis_params['npol'] = corrparams.npol
    vis_params['nfint'] = corrparams.nfint
    vis_params['ntint'] = corrparams.ntint
    return vis_params


def parse_visibility_parameters(
        params: dict, tstart: 'astropy.time.Time', ntint: int) -> dict:
    """Parse details about the data to be extracted.

    Parameters
    ----------
    params : dict
        Parameters passed by the user definining the correlator and T3 system.
    tstart : astropy.time.Time
        Start time of the correlated visibility file.
    ntint : int
        The number of time bins integrated together in the correlated data
        compared to the native resolution.

    Returns
    -------
    dict
        Parameters that describe the visibilities to be read in.
    """
    antenna_order = params['antennas']
    fobs = params['f0_GHz']+params['deltaf_MHz']*1e-3*np.arange(params['nchan'])
    nant = len(antenna_order)
    nbls = (nant*(nant+1))//2

    # Visibilities have already been integrated in time, so account for that when
    # determining the observation times.
    tsamp = params['deltat_s']*ntint*u.s
    tobs = tstart + (np.arange(params['nsubint']//ntint)+0.5)*tsamp
    print("Time...",tobs.isot)
    
    # Get baselines
    blen, bname = get_blen(antenna_order)
    # Get indices for baselines to the reference antenna
    refidxs = []
    refant = str(antenna_order[0])
    for i, bn in enumerate(bname):
        if refant in bn.split('-'):
            refidxs += [i]

    print(params['outrigger_delays'])
    cable_delays = get_cable_delays(params['outrigger_delays'], bname)

    vis_params = {
        # Baselines
        'antenna_order': antenna_order,
        'blen': blen,
        'bname': bname,
        'nbls': nbls,
        'refidxs': refidxs,
        'antenna_cable_delays': params['outrigger_delays'],
        'baseline_cable_delays': cable_delays,
        'snapdelays': params['snapdelays'],
        'ant_itrf': params['ant_itrf'],
        'nants_telescope': params['nants_telescope'],
        # Time
        'tsamp': tsamp,
        'tobs': tobs,
        # Frequency
        'fobs': fobs,
        'corr_ch0': params['ch0'],
        'nchan_corr': params['nchan_corr']}
    return vis_params


def get_cable_delays(outrigger_delays: dict, bname: list) -> np.ndarray:
    """Calculate cable delays from the measured outrigger cable delays.

    Antenna names in both input parameters are indexed at 1.

    Parameters
    ----------
    outrigger_delays : dict
        The delay for each antenna.  Missing keys have 0 delay.
    bname : list
        The name of each baseline.

    Returns
    -------
    np.ndarray
        The cable delay for each baseline in bname.
    """
    print(outrigger_delays.keys())
    delays = np.zeros(len(bname), dtype=np.int)
    for i, bn in enumerate(bname):
        ant1, ant2 = bn.split('-')
        delays[i] = outrigger_delays.get(int(ant1), 0)-\
                    outrigger_delays.get(int(ant2), 0)
    return delays


def get_blen(antennas: list) -> tuple:
    """Gets the baseline lengths for a subset of antennas.

    Parameters
    ----------
    antennas : list
        The antennas used in the array.

    Returns
    -------
    blen : array
        The ITRF coordinates of all of the baselines.
    bname : list
        The names of all of the baselines.
    """
    ant_itrf = get_itrf(
        latlon_center=(ct.OVRO_LAT*u.rad, ct.OVRO_LON*u.rad, ct.OVRO_ALT*u.m)
    ).loc[antennas]
    xx = np.array(ant_itrf['dx_m'])
    yy = np.array(ant_itrf['dy_m'])
    zz = np.array(ant_itrf['dz_m'])
    # Get uvw coordinates
    nants = len(antennas)
    nbls = (nants*(nants+1))//2
    blen = np.zeros((nbls, 3))
    bname = []
    k = 0
    for i in range(nants):
        for j in range(i, nants):
            blen[k, :] = np.array([
                xx[i]-xx[j],
                yy[i]-yy[j],
                zz[i]-zz[j]
            ])
            bname += ['{0}-{1}'.format(
                antennas[i],
                antennas[j]
            )]
            k += 1
    return blen, bname


class ProtectedVariable():
    """A protected variable with variable and lock managed by dask"""

    def __init__(self, name, initial_value=None):
        self.name = name
        self.variable = Variable(name)
        self.lock = Lock(name)
        self.variable.set(initial_value)

    @contextmanager
    def get_lock(self):
        try:
            assert self.lock.acquire()
            yield self.variable
        finally:
            self.lock.release()

    def get(self):
        return self.variable.get()

    def set(self, value):
        self.variable.set(value)
