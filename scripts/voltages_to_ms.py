"""
Convert voltage files to measurement sets.
"""
import os
from multiprocessing import Process, Manager
import argparse
from dsaT3.uvh5_to_ms import uvh5_to_ms
from dsaT3.voltages_to_ms import *

def voltages_to_ms(candname: str, datestring: str, ntint: int, start_offset: int, end_offset: int,
                   dispersion_measure: float=None, full_pol: bool=False, continuum_source: bool=False) -> None:
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
    dispersion_measure: float
        Overrides the dispersion measure in the json header file.
    full_pol: bool
        If True, no frequency averaging is done, and 4 polarizations are written out.
        If False, frequency averaging x 8 is done, and XX and YY are written out.
    """

    print(f'Continuum Source: {continuum_source}')


    start_offset, end_offset = set_default_if_unset(start_offset, end_offset)
    system_setup = initialize_system()
    cand = initialize_candidate(candname, datestring, system_setup, dispersion_measure)
    corrparams = initialize_correlator(full_pol, ntint, cand, system_setup)
    uvh5params = initialize_uvh5(corrparams, cand, system_setup)

    # Initialize the process manager, locks, values, and queues
    manager = Manager()
    ncorrfiles = manager.Value('i', lock=True)
    ncorrfiles.value = 0
    declination = manager.Value(float, lock=True)
    declination.value = None
    rsync_queue = manager.Queue()
    corr_queue = manager.Queue()
    uvh5_queue = manager.Queue()

    # Generate the table needed by the correlator
    get_declination_etcd = process_join(generate_declination_component(
        declination, cand.time))
    _ = get_declination_etcd()

    generate_delay_table(uvh5params.visparams, corrparams.reftime, declination.value)

    rsync_all_files = pipeline_component(
        generate_rsync_component(cand.local),
        rsync_queue,
        corr_queue)

    correlate = pipeline_component(
        generate_correlate_component(
            cand.dm, corrparams.ntint, system_setup.corr_ch0_MHz,
            corrparams.npol, ncorrfiles),
        corr_queue,
        uvh5_queue)

    write_uvh5 = pipeline_component(
        generate_uvh5_component(
            cand.name, system_setup.corrdir, declination, uvh5params.visparams,
            start_offset, end_offset, ncorrfiles),
        uvh5_queue)

    processes = [
        Process(target=rsync_all_files, daemon=True),
        Process(target=correlate, daemon=True),
        Process(target=write_uvh5, daemon=True)]

    # Populate rsync queue (the entry queue in our pipeline)
    # We only need to do this if the corresponding hdf5 file hasn't
    # been created already.
    for i, filename in enumerate(cand.voltagefiles):
        if not os.path.exists(uvh5params.files[i]):
            rsync_queue.put([filename, corrparams.files[i]])
    rsync_queue.put('END')

    # Start all processes in the pipeline
    for proc in processes:
        proc.start()

    # Wait for all processes to finish
    for proc in processes:
        proc.join()

    # Convert uvh5 files to a measurement set
    msname = f'{system_setup.msdir}{candname}'
    ntbins = None if continuum_source else 8
    uvh5_to_ms(cand.name, cand.time, uvh5params.files, msname, corrparams.reftime, ntbins)

    # Remove hdf5 files from disk
    # for hdf5file in uvh5params.files:
    #     os.remove(hdf5file)

def set_default_if_unset(start_offset: int, end_offset: int) -> tuple:
    """Set default `start_offset` and `end_offset` if they aren't set."""
    if start_offset < 0:
        start_offset = None
    if end_offset < 0:
        end_offset = None
    return start_offset, end_offset

def parse_commandline_arguments() -> "argparse.Namespace":
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description='Correlate candidate voltage files.')
    parser.add_argument(
        'candname',
        type=str,
        help='unique candidate name')
    parser.add_argument(
        '--datestring',
        type=str,
        help='datestring of archived candidate',
        nargs='?',
        default='current')
    parser.add_argument(
        '--ntint',
        type=int,
        nargs='?',
        default=8,
        help='number of native time bins to integrate during correlation')
    parser.add_argument(
        '--startoffset',
        type=int,
        nargs='?',
        default=1716,
        help='nbins from beginning of correlated data to start writing to ms')
    parser.add_argument(
        '--stopoffset',
        type=int,
        nargs='?',
        default=4788,
        help='number of bins from end of correlation to write to ms')
    parser.add_argument(
        '--dm',
        type=float,
        nargs='?',
        help='dispersion measure to use instead of value in header file')
    parser.add_argument('--continuum', action='store_true')
    parser.set_defaults(continuum=False)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    ARGS = parse_commandline_arguments()
    voltages_to_ms(ARGS.candname, ARGS.datestring, ntint=ARGS.ntint,
                   start_offset=ARGS.startoffset, end_offset=ARGS.stopoffset,
                   dispersion_measure=ARGS.dm, continuum_source=ARGS.continuum)
