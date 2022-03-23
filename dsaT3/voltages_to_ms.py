"""Utilties for converting T3 voltage dumps to measurement sets."""

from collections import namedtuple
import time
import re
import os
from functools import wraps
import subprocess
from multiprocessing import Process
import queue
from pkg_resources import resource_filename
import numpy as np
import astropy.units as u
from antpos.utils import get_itrf
from dsautils.coordinates import get_declination, get_elevation
import dsacalib.constants as ct
from dsaT3.generate_uvh5 import calculate_uvw_and_geodelay, get_total_delay
from dsaT3.utils import rsync_file, load_params, get_tstart_from_json, get_DM_from_json
from dsaT3.generate_uvh5 import generate_uvh5

__all__ = ['pipeline_component', 'generate_rsync_component',
           'generate_correlate_component', 'generate_uvh5_component',
           'process_join', 'generate_declination_component', 'generate_delay_table',
           'initialize_system', 'initialize_candidate', 'initialize_correlator',
           'initialize_uvh5', 'initialize_vis_params', 'parse_visibility_parameters',
           'get_cable_delays', 'get_blen']

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')

def pipeline_component(targetfn, inqueue, outqueue=None):
    """Generate a component of the pipeline."""
    @wraps(targetfn)
    def inner():
        """Process data from a queue."""
        done = False
        while not done:
            try:
                item = inqueue.get()
                assert item
            except (queue.Empty, AssertionError) as _:
                time.sleep(10)
                continue
            except (EOFError, BrokenPipeError) as exc:
                print(f'Error when accessing inqueue in {targetfn.__name__}')
                done = True
            if item == 'END':
                done = True
            else:
                item = targetfn(item)
            try:
                if outqueue is not None:
                    outqueue.put(item)
            except (EOFError, BrokenPipeError) as exc:
                print(f'Error when accessing outqueue in {targetfn.__name__}')
                done = True
    return inner

def generate_rsync_component(local: bool) -> "Callable":
    """Generate an rsync function."""

    def rsync_all_files(item):
        """Rsync or copy a file."""
        srcfile, vfile = item
        if not os.path.exists(vfile):
            if not local:
                rsync_file(srcfile, vfile)
            else:
                os.symlink(srcfile, vfile)
        return vfile
    return rsync_all_files

def generate_correlate_component(
        dispersion_measure: float, ntint: int, corr_ch0: dict, npol: int,
        ncorrfiles: "Manager().Value") -> "Callable":
    """Generate a correlator function."""

    def correlate(vfile):
        """Correlate a file."""
        while ncorrfiles.value > 2:
            time.sleep(10)

        corr = re.findall('corr\d\d', vfile)[0]
        if not os.path.exists('{0}.corr'.format(vfile)):
            first_channel_MHz = corr_ch0[corr]
            command = (
                '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/toolkit_dev '
                f'-i {vfile} -o {vfile}.corr -t {ntint} -c {first_channel_MHz} '
                f'-d delays.dat {"" if npol==4 else "-a"}')
            if dispersion_measure is not None:
                command += f' -m {dispersion_measure}'
            print(command)
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                shell=True)
            proc_stdout = str(process.communicate()[0].strip())
            print(proc_stdout)

        corrfile = f'{vfile}.corr'

        with ncorrfiles.get_lock():
            ncorrfiles.value += 1


        return corrfile

    return correlate

def generate_uvh5_component(
        candname: str, corrdir: str, declination: "Manager().Value", vis_params: dict,
        start_offset: int, end_offset: int, ncorrfiles: "Manager().Value") -> "Callable":
    """Generate a uvh5 writer."""

    def write_uvh5(corrfile):
        """Write correlated data to a uvh5 file."""
        print(corrfile, type(corrfile))
        _uvh5name = generate_uvh5(
            f'{corrdir}/{candname}',
            declination.value*u.deg,
            corrfile=corrfile,
            vis_params=vis_params,
            start_offset=start_offset,
            end_offset=end_offset
        )

        os.remove(corrfile)
        try:
            with ncorrfiles.get_lock():
                ncorrfiles.value -= 1
        except (EOFError, BrokenPipeError) as exc:
            print('Error when writing to ncorrfiles frome write_uvh5')

    return write_uvh5

def process_join(targetfn):
    """Generate a function that starts and process and waits for it to complete."""

    def inner():
        """Start a process and then wait for it to complete."""
        process = Process(
            target=targetfn,
            daemon=True)
        process.start()
        process.join()
        return process
    return inner

def generate_declination_component(
        declination: "Value", tstart: "astropy.time.Time") -> "Callable":
    """Generate a pipeline component to get the declination from etcd."""

    def get_declination_etcd():
        """Look up the declination from etcd."""
        declination.value = get_declination(
            get_elevation(tstart)
        ).to_value(u.deg)

    return get_declination_etcd

def generate_delay_table(vis_params, reftime, declination):
    """Generate a table of geometric and cable delays for the correlator."""
    _buvw, ant_bw = calculate_uvw_and_geodelay(vis_params, reftime, declination*u.deg)
    total_delay = get_total_delay(
        vis_params['baseline_cable_delays'], ant_bw, vis_params['bname'],
        vis_params['antenna_order'])
    total_delay_string = '\n'.join(total_delay.flatten().astype('str'))+'\n'
    with open("delays.dat", "w", encoding='utf-8') as f:
        f.write(total_delay_string)

def initialize_system():
    """Set system parameters needed in voltage to ms converter."""
    params = load_params(PARAMFILE)
    SystemSetup = namedtuple(
        "SystemSetup", "T3dir archivedir corrdir msdir start_time_offset reffreq_GHz corr_ch0_MHz")
    start_time_offset = params['burst_start_s']*u.s
    corr_ch0_MHz = {corr: 1e3*params['f0_GHz']+params['deltaf_MHz']*corr_ch0
                    for corr, corr_ch0 in params['ch0'].items()}

    system_setup = SystemSetup(
        params['T3dir'], params['archivedir'], params['corrdir'], params['msdir'],
        start_time_offset, params['reffreq_GHz'], corr_ch0_MHz)
    return system_setup

def initialize_candidate(candname, datestring, system_setup, dispersion_measure=None):
    """Set candidate parameters using information in the json header."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())

    if datestring == 'current':
        local = False
        headerfile = f'{system_setup.t3dir}/{candname}.json'
        voltagefiles = [f'{corr}.sas.pvt:/home/ubuntu/data/{candname}_data.out'
                        for corr in corrlist]

    else:
        local = True
        headerfile = f'{system_setup.archivedir}/{datestring}/{candname}.json'
        voltagefiles = [f'{system_setup.archivedir}/{datestring}/{corr}_{candname}_data.out'
                        for corr in corrlist]

    tstart = get_tstart_from_json(headerfile)

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
    reftime = cand.time + system_setup.start_time_offset
    npol = 4 if fullpol else 2
    nfint = 1 if fullpol else 8
    corrfiles = [f'{system_setup.corrdir}/{corr}_{cand.name}_data.out' for corr in corrlist]

    CorrelatorParameters = namedtuple('Correlator', 'reftime npol nfint ntint files')
    correlator_params = CorrelatorParameters(reftime, npol, nfint, ntint, corrfiles)
    return correlator_params

def initialize_uvh5(corrparams, cand, system_setup):
    """Set parameters for writing uvh5 files."""
    corrlist = list(system_setup.corr_ch0_MHz.keys())
    uvh5files = [f'{system_setup.corrdir}/{cand.name}_{corr}.hdf5' for corr in corrlist]
    visparams = initialize_vis_params(corrparams, cand)

    UVH5Parameters = namedtuple('UVH5', 'files visparams')
    uvh5_params = UVH5Parameters(uvh5files, visparams)
    return uvh5_params

def initialize_vis_params(corrparams, cand):
    """Set parameters that describe the visbilities produced by the correlator."""
    T3params = load_params(PARAMFILE)
    vis_params = parse_visibility_parameters(T3params, cand.time, corrparams.ntint)
    vis_params['tref'] = corrparams.reftime
    vis_params['npol'] = corrparams.npol
    vis_params['nfint'] = corrparams.nfint
    vis_params['ntint'] = corrparams.ntint
    return vis_params

def parse_visibility_parameters(params: dict, tstart: 'astropy.time.Time', ntint: int) -> dict:
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

    # Get baselines
    blen, bname = get_blen(antenna_order)
    # Get indices for baselines to the reference antenna
    refidxs = []
    refant = str(antenna_order[0])
    for i, bn in enumerate(bname):
        if refant in bn.split('-'):
            refidxs += [i]

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
