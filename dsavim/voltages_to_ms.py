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
from dsavim.generate_uvh5 import calculate_uvw_and_geodelay, get_total_delay, calculate_uvw_and_geodelay_radec
from dsavim.utils import rsync_file, get_tstart_from_json, get_DM_from_json, load_params
from dsavim.generate_uvh5 import generate_uvh5

from dsamfs.fringestopping import generate_fringestopping_table

__all__ = [
    'pipeline_component', 'rsync_component', 'correlate_component',
    'uvh5_component',
    'generate_delay_table', 'generate_delay_table_radec', 'initialize_system', 'initialize_candidate', 'initialize_correlator',
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
    _corrlist = {'sb00':'lxd110h03','sb01':'lxd110h04','sb02':'lxd110h05','sb03':'lxd110h06','sb04':'lxd110h07','sb05':'lxd110h08','sb06':'lxd110h10','sb07':'lxd110h11','sb08':'lxd110h12','sb09':'lxd110h14','sb10':'lxd110h15','sb11':'lxd110h16','sb12':'lxd110h18','sb13':'lxd110h19','sb14':'lxd110h21','sb15':'lxd110h22'}
    corr = _corrlist[sb]
    
    if not os.path.exists(f"{vfile}.corr"):
        first_channel_MHz = corr_ch0[corr]
        # dsa110-bbproc `toolkit` visibilities mode (2026-07: replaces
        # the retired toolkit_dev binary; same -i/-o/-t/-c/-d/-a/-m CLI
        # and float32 output, but reads the TRUE fragment length
        # instead of a hardcoded 4.5 GB, and -m dedisperses host-side).
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

def generate_delay_table_radec(vis_params, reftime, rajd, declination):
    """Generate a table of geometric and cable delays for the correlator, with RA and DEC."""

    _buvw, ant_bw = calculate_uvw_and_geodelay_radec(vis_params, reftime, (rajd*np.pi/180.)*u.rad, (declination*np.pi/180.)*u.rad)
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


def initialize_candidate(candname, system_setup, dispersion_measure=None,
                         time_shift_s=0.0):
    """Set candidate parameters from the dsa110-rt (M8) archive.

    2026-07 rewrite: the T3-era ``T2_<cand>.json`` is gone; burst time,
    DM and the voltage-file inventory come from
    :func:`dsavim.metadata.load_burst` (Level3 c2 record + per-subband
    staging manifests). ``cand.time`` is the FILE START (exact, from
    the dump manifests) and ``cand.nsamples`` the per-fragment native
    sample count — the T3 code hardcoded both (burst_start_s /
    nsubint yaml constants)."""
    from dsavim.metadata import load_burst

    burst = load_burst(candname, dm=dispersion_measure,
                       time_shift_s=time_shift_s)
    Candidate = namedtuple(
        'Candidate',
        'name time dm voltagefiles headerfile local burst nsamples')
    return Candidate(
        candname, burst.file_start, burst.dm,
        [str(f) for f in burst.voltage_files],
        None, True, burst, burst.n_native_samples)


def initialize_correlator(fullpol, ntint, cand, system_setup):
    """Set correlator parameters.

    ``reftime`` (the phase/fringestop centre time downstream) is the
    burst arrival — from the M8 manifests via ``cand.burst`` (the T3
    code used file_start + a fixed 0.4999 s)."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    sblist = ['sb00','sb01','sb02','sb03','sb04','sb05','sb06','sb07',
              'sb08','sb09','sb10','sb11','sb12','sb13','sb14','sb15']
    import astropy.units as u
    reftime = cand.time + cand.burst.burst_offset_s * u.s
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

    # M8: per-fragment sample count comes from the dump manifests, not
    # the yaml nsubint constant (T3 files were always 61440 samples;
    # M8 fragments are n_blocks_written x 4096, typically 94208).
    T3params = dict(T3params)
    T3params['nsubint'] = cand.nsamples
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
    delays = np.zeros(len(bname), dtype=int)
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
