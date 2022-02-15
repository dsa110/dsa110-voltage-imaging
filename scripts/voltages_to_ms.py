"""
Convert voltage files to measurement sets.
"""
from types import MappingProxyType
import re
import os
import glob
import subprocess
from multiprocessing import Process, Manager
import multiprocessing
import queue
import argparse
import time
from pkg_resources import resource_filename
import astropy.units as u
from dsautils.coordinates import get_declination, get_elevation
from dsaT3.uvh5_to_ms import uvh5_to_ms
from dsaT3.utils import rsync_file, load_params, get_tstart_from_json, get_DM_from_json
from dsaT3.generate_uvh5 import generate_uvh5, parse_visibility_parameters

NPROC = 2
PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
T3PARAMS = load_params(PARAMFILE)
CORR_LIST = list(T3PARAMS['ch0'].keys())
GEN_DELAY_SCRIPT = ('/home/ubuntu/anaconda3/envs/dana/bin/python '
                    '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/gen_delays.py')
BURST_START_S = 1907*262.144e-6

def voltages_to_ms(candname: str, datestring: str, ntint: int, start_offset: int, end_offset: int,
                   full_pol: bool=False) -> None:
    """
    Correlate voltage files and convert to a measurement set.

    Parameters
    ----------
    candname : str
        The unique name of the candidate.
    datestring : str
        The datestring the observation is archived under. Use 'current' if the
        data is from the current, unarchived observing run.
    ntint : int
        The number of time samples to integrate together during correlation.
    nfint : int
        The number of frequency channels to integrate together after removing
        outrigger delays.
    start_offset : int
        The number of time samples (after correlation) to offset the start of
        the measurement set by.  If not provided, the entire time is converted
        to a measurement set.
    end_offset : int
        The last time sample (after correlation) to write to the measurement
        set. If not provide,d the entire time is converted to a measurement
        set.
    """

    start_offset, end_offset = set_default_if_unset(start_offset, end_offset)
    filenames, headername, rsync = get_input_file_locations(candname, datestring)
    outnames, hdf5names = get_output_file_locations(candname)

    # Get additional parameters that describe the correlation
    tstart = get_tstart_from_json(headername)
    dispersion_measure = get_DM_from_json(headername)
    if dispersion_measure < 1:
        dispersion_measure = None

    vis_params = parse_visibility_parameters(T3PARAMS, tstart, ntint)
    vis_params['tref'] = tstart+BURST_START_S
    vis_params['npol'] = 4 if full_pol else 2
    vis_params['nfint'] = 1 if full_pol else 8
    vis_params['ntint'] = ntint
    corr_ch0_MHz = {key: 1e3*vis_params['fobs'][value] for key, value in vis_params['corr_ch0'].items()}
    corr_ch0_MHz_safe = MappingProxyType(corr_ch0_MHz) # Thread-safe mapping
    vis_params_safe = MappingProxyType(vis_params)

    # Initialize the process manager, locks, values, and queues
    manager = Manager()
    ncorrfiles = manager.Value('i', 0)
    ncorrfiles_lock = manager.Lock()
    declination = manager.Value(float, None)
    declination_lock = manager.Lock()
    rsync_queue = manager.Queue()
    corr_queue = manager.Queue()
    uvh5_queue = manager.Queue()

    # Generate the table needed by the correlator
    process = Process(
        target=get_declination_handler,
        args=(declination, declination_lock, tstart),
        daemon=True)
    process.start()
    process.join()

    # TODO: Add this to generate_uvh5 and import from there
    # TODO: Generate this when correlating
    command = f'{GEN_DELAY_SCRIPT} {headername} {declination.value}'
    process = subprocess.Popen(
        command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        shell=True)
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)

    # Populate rsync queue (the entry queue in our pipeline)
    # We only need to do this if the corresponding hdf5 file hasn't
    # been created already.
    for i, filename in enumerate(filenames):
        if not os.path.exists(hdf5names[i]):
            rsync_queue.put([filename, outnames[i]])
    rsync_queue.put('END')

    # Define the processes that make up the pipeline
    processes = []
    # Rsync processes copy data files to the local directory for processing
    processes += [Process(
        target=rsync_handler,
        args=(rsync_queue, corr_queue, rsync),
        daemon=True)]
    # Corr processes correlate the data
    processes += [Process(
        target=corr_handler,
        args=(ntint, corr_ch0_MHz_safe, full_pol, corr_queue,
              uvh5_queue, ncorrfiles, ncorrfiles_lock),
        daemon=True)]
    for i in range(NPROC):
        # UVH5 processes convert the correlated data to uvh5 format
        processes += [Process(
            target=uvh5_handler,
            args=(candname, declination, vis_params_safe,
                  start_offset, end_offset, uvh5_queue, ncorrfiles, ncorrfiles_lock),
            daemon=True)]

    # Start all processes in the pipeline
    for proc in processes:
        proc.start()

    # Wait for rsync process to finish
    processes[0].join()
    print('rsync process done')
    corr_queue.put('END')

    # Wait for corr process to finish
    processes[1].join()
    print('corr process done')
    for i in range(NPROC):
        uvh5_queue.put('END')

    # Wait for uvh5 processes to finish
    for proc in processes[2:]:
        proc.join()
        print('A uvh5 process joined.')
    print('All uvh5 processes done')

    # Convert uvh5 files to a measurement set
    # This can also be done with a queue
    hdf5files = sorted(glob.glob(f'{T3PARAMS["corrdir"]}/{candname}_corr??.hdf5'))
    uvh5_to_ms(candname, tstart, dispersion_measure, hdf5files, f'{T3PARAMS["msdir"]}/{candname}')

    # Remove hdf5 files from disk
    for hdf5file in hdf5files:
        os.remove(hdf5file)

def rsync_handler(rsync_queue: "Manager().Queue", corr_queue: "Manager().Queue",
                  rsync: bool) -> None:
    """Rsyncs files, then updates `corr_queue`."""
    rsync_done = False
    while not rsync_done:
        try:
            item = rsync_queue.get()
        except queue.Empty:
            time.sleep(10)
        else:
            if item == 'END':
                rsync_done = True
                continue
            srcfile, vfile = item
            if not os.path.exists(vfile):
                if rsync:
                    rsync_file(srcfile, vfile)
                else:
                    os.symlink(srcfile, vfile)
            corr_queue.put(vfile)

def corr_handler(ntint: int, corr_ch0: dict, full_pol: bool, corr_queue: "Manager().Queue",
                 uvh5_queue: "Manager().Queue", ncorrfiles: "Manager().Value",
                 ncorrfiles_lock: "Manager().Lock") -> None:
    """Correlates data using T3 cpu correlator."""
    corr_done = False
    while not corr_done:
        if ncorrfiles.value > NPROC-1:
            time.sleep(10)
            continue
        try:
            vfile = corr_queue.get()
        except queue.Empty:
            time.sleep(10)
        else:
            if vfile == 'END':
                corr_done = True
                continue
            with ncorrfiles_lock:
                ncorrfiles.value += 1
            corr = re.findall('corr\d\d', vfile)[0]
            if not os.path.exists('{0}.corr'.format(vfile)):
                first_channel_MHz = corr_ch0[corr]
                command = (
                    '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/toolkit '
                    f'-i {vfile} -o {vfile}.corr -t {ntint} -c {first_channel_MHz} '
                    f'-d delays.dat {"" if full_pol else "-a"}')
                print(command)
                process = subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    shell=True)
                proc_stdout = str(process.communicate()[0].strip())
                print(proc_stdout)

            corrfile = '{0}.corr'.format(vfile)
            uvh5_queue.put(corrfile)

def get_declination_handler(declination: "Manager().Value", declination_lock: "Manager().Lock",
                    tstart: "astropy.time.Time"):
    with declination_lock:
        if declination.value is None:
            declination.value = get_declination(
                get_elevation(tstart)
            ).to_value(u.deg)

def uvh5_handler(candname: str, declination: "Manager().Value",
                 vis_params: dict, start_offset: int,
                 end_offset: int, uvh5_queue: "Manager().Queue", ncorrfiles: "Manager().Value",
                 ncorrfiles_lock: "Manager().Lock") -> None:
    """Convert correlated data to uvh5."""
    proc = multiprocessing.current_process()
    uvh5_done = False

    while not uvh5_done:
        try:
            corrfile = uvh5_queue.get()
        except queue.Empty:
            time.sleep(10)
        else:
            if corrfile == 'END':
                print('proc {0} is setting uvh5_done'.format(proc.pid))
                uvh5_done = True
                continue
            uvh5name = generate_uvh5(
                '{0}/{1}'.format(T3PARAMS['corrdir'], candname),
                declination.value*u.deg,
                corrfile=corrfile,
                vis_params=vis_params,
                start_offset=start_offset,
                end_offset=end_offset
            )
            print(uvh5name)
            os.remove(corrfile)
            with ncorrfiles_lock:
                ncorrfiles.value -= 1
    print('{0} exiting'.format(proc.pid))

def set_default_if_unset(start_offset: int, end_offset: int) -> tuple:
    """Set default `start_offset` and `end_offset` if they aren't set."""
    if start_offset < 0:
        start_offset = None
    if end_offset < 0:
        end_offset = None
    return start_offset, end_offset

def get_input_file_locations(candname: str, datestring: str) -> tuple:
    """Get the locations of input data and header files."""
    if datestring == 'current':
        rsync = True
        filenames = [
            f'{corr}.sas.pvt:/home/ubuntu/data/{candname}_data.out'
            for corr in CORR_LIST]
        headername = f'{T3PARAMS["T3dir"]}/{candname}.json'
    else:
        rsync = False
        filenames = [
            f'{T3PARAMS["archivedir"]}/{datestring}/{corr}_{candname}_data.out'
            for corr in CORR_LIST]
        headername = f'{T3PARAMS["archivedir"]}/{datestring}/{candname}.json'
    return filenames, headername, rsync

def get_output_file_locations(candname: str) -> tuple:
    """Determine where to put intermediate output files (correlated and uvh5)"""
    outnames = [f'{T3PARAMS["corrdir"]}/{corr}_{candname}_data.out' for corr in CORR_LIST]
    hdf5names = [f'{T3PARAMS["corrdir"]}/{candname}_{corr}.hdf5' for corr in CORR_LIST]
    return outnames, hdf5names

def parse_commandline_arguments() -> "argparse.Namespace":
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Correlate candidate voltage files.'
    )
    parser.add_argument(
        'candname',
        type=str,
        help='unique candidate name'
    )
    parser.add_argument(
        '--datestring',
        type=str,
        help='datestring of archived candidate',
        nargs='?',
        default='current'
    )
    parser.add_argument(
        '--ntint',
        type=int,
        nargs='?',
        default=8,
        help='number of native time bins to integrate during correlation'
    )
    parser.add_argument(
        '--startoffset',
        type=int,
        nargs='?',
        default=1716,
        help='nbins from beginning of correlated data to start writing to ms'
    )
    parser.add_argument(
        '--stopoffset',
        type=int,
        nargs='?',
        default=4788, #3252, #2484,
        help='number of bins from end of correlation to write to ms'
    )
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    ARGS = parse_commandline_arguments()
    voltages_to_ms(ARGS.candname, ARGS.datestring, ntint=ARGS.ntint,
                   start_offset=ARGS.startoffset, end_offset=ARGS.stopoffset)
