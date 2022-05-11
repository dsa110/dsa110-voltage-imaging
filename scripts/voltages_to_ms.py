"""
Convert voltage files to measurement sets.
"""
print("Importing stuff")
import os
import argparse
from functools import partial

import dask
dask.config.set(scheduler="synchronous")
from dask.distributed import Client, Variable, wait

from astropy.time import Time

import dsautils.cnf as dsc

from dsavim.uvh5_to_ms import uvh5_to_ms
from dsavim.voltages_to_ms import *


def voltages_to_ms(
        candname: str, datestring: str, ntint: int, start_offset: int, end_offset: int,
        dispersion_measure: float = None, full_pol: bool = False) -> None:
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

    # The following needs to happen in a subprocess
    system_setup = initialize_system()

    cand = initialize_candidate(candname, datestring, system_setup, dispersion_measure)
    corrparams = initialize_correlator(full_pol, ntint, cand, system_setup)

    # TODO: Look up outrigger delays in influx instead of using these hardcoded values
    if cand.time < Time("2022-03-17"):
        outrigger_delays = {
            100:  1256,
            101:  1106,
            102:  1054,
            103:  -212,
            104:  2129,
            105: 11992,
            106:  8676,
            107:  9498,
            108: 10310,
            109: 11438,
            110: 15450,
            111: 14414,
            112: 15158,
            113: 16622,
            114: 18742,
            115: 20298,
            116:  3676,
            117:  4376}
    else:
        outrigger_delays = None

    uvh5params = initialize_uvh5(corrparams, cand, system_setup, outrigger_delays)

    reftime = get_reftime()
    declination = get_declination_etcd(cand.time)

    corr2uvh5 = UVH5Generator(candname=cand.name, corrdir=system_setup.corrdir,
            declination=declination, vis_params=uvh5params.visparams, start_offset=start_offset,
            end_offset=end_offset)

    ant_bw = uvh5_gen.calculate_uvw_and_geodelay(tobs=reftime)
    generate_delay_table(ant_bw, corr2uvh5.visparams)

    client = Client(name="dsavim", n_workers=2, scheduler_port=9999, dashboard_address="localhost:9998")

    rsync_all_files = partial(rsync_component, local=cand.local)
    correlate = partial(
        correlate_component, dispersion_measure=cand.dm, ntint=corrparams.ntint,
        corr_ch0=system_setup.corr_ch0_MHz, npol=corrparams.npol)

    futures = []
    last_corr_future = ""
    for i, filename in enumerate(cand.voltagefiles):
        if not os.path.exists(uvh5params.files[i]):
            rsync_future = client.submit(rsync_all_files, [filename, corrparams.files[i]])
            # TODO: Have one single worker with gpu resources to process correlate instead of
            # making it depend on the previous process finishing
            corr_future = client.submit(correlate, rsync_future, last_corr_future)
            corr2uvh5_future = client.submit(corr2uvh5.process, corr_future)
            futures += [rsync_future, corr_future, corr2uvh5_future]
            last_corr_future = corr_future
    wait(futures)
    client.close()

    # Convert uvh5 files to a measurement set
    msname = f"{system_setup.msdir}{candname}"
    uvh5_to_ms(
        cand.name, cand.time, uvh5params.files, msname, corrparams.reftime,
        template_path=None)

    # Remove hdf5 files from disk
    for hdf5file in uvh5params.files:
        os.remove(hdf5file)

def set_default_if_unset(start_offset: int, end_offset: int) -> tuple:
    """Set default `start_offset` and `end_offset` if they aren't set."""
    if start_offset < 0:
        start_offset = None
    if end_offset < 0:
        end_offset = None
    return start_offset, end_offset

def get_reftime() -> "Time":
    """Get the reference time used in the real-time correlator."""
    conf = dsc.Conf()
    refmjd = conf.get('fringe')['refmjd']
    reftime = Time(refmjd, format='mjd')
    return reftime

def parse_commandline_arguments() -> 'argparse.Namespace':
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description="Correlate candidate voltage files.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        'candname',
        type=str,
        help="unique candidate name")
    parser.add_argument(
        '--datestring',
        type=str,
        help="datestring of archived candidate",
        nargs='?',
        default='current')
    parser.add_argument(
        '--ntint',
        type=int,
        nargs='?',
        default=32,
        help="number of native time bins to integrate during correlation")
    parser.add_argument(
        '--startoffset',
        type=int,
        nargs='?',
        default=444,
        help="nbins from beginning of correlated data to start writing to ms")
    parser.add_argument(
        '--stopoffset',
        type=int,
        nargs='?',
        default=508,
        help='number of bins from end of correlation to write to ms')
    parser.add_argument(
        '--dm',
        type=float,
        nargs='?',
        help="dispersion measure to use instead of value in header file")

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    print("Running main program")
    ARGS = parse_commandline_arguments()
    voltages_to_ms(
        ARGS.candname, ARGS.datestring, ntint=ARGS.ntint, start_offset=ARGS.startoffset,
        end_offset=ARGS.stopoffset, dispersion_measure=ARGS.dm)
