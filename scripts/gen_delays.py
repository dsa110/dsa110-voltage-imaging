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
from astropy.coordinates import SkyCoord
import astropy.units as u

import dsautils.cnf as dsc

from dsavim.uvh5_to_ms import uvh5_to_ms
from dsavim.voltages_to_ms import *
import sys
import numpy as np

def voltages_to_ms(
        candname: str, ntint: int, start_offset: int, end_offset: int,
        dispersion_measure: float = None, full_pol: bool = False, weights: str = None, obsra: float = None, obsdec: float = None, burstra: float = None, burstdec: float = None) -> None:
    """
    Correlate voltage files and convert to a measurement set.

    Parameters
    ----------
    candname : str
        The unique name of the candidate.
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
    weights : str
        (optional) the datestring of any weights to apply pre-correlation
    """

    start_offset, end_offset = set_default_if_unset(start_offset, end_offset)

    # The following needs to happen in a subprocess
    system_setup = initialize_system()

    cand = initialize_candidate(candname, system_setup, dispersion_measure)
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
    rajd = obsra
    declination = obsdec#get_declination_etcd(cand.time)

    f = open("baselines.txt","w")
    for i in np.arange(len(uvh5params.visparams['bname'])):
        f.write(f"{i} {uvh5params.visparams['bname'][i]}\n")
    f.close()
    #sys.exit(1)
        

    
    #generate_delay_table(uvh5params.visparams, reftime, declination, use_fs=False)
    #os.system("mv delays.dat delays_meridian.dat")

    # now generate antenna delay table for offset to burst
    #phasecenter = SkyCoord(obsra,obsdec,unit=(u.deg,u.deg))
    #burst = SkyCoord(burstra,burstdec,unit=(u.deg,u.deg))
    #record_burst_delays(phasecenter,burst)
    
    generate_delay_table_radec(uvh5params.visparams, cand.time, rajd, declination)
    os.system("mv delays.dat delays_meridian.dat")
    
def record_burst_delays(phasecenter,burst):

    dra, ddec = phasecenter.spherical_offsets_to(burst)
    print(dra.value,ddec.value) # these are in deg
    drar = dra.value*np.pi/180.
    ddecr = ddec.value*np.pi/180.
    

    # RA and DEC are -offset*antpos/c. Then in toolkit it's 2*pi*VAL*freqs as the phasor
    _antpos = np.load("antpos.npz")
    antpos = _antpos["antpos"] # (2,96) array
    delays = np.zeros(96)
    for i in np.arange(96):
        delays[i] = 1e9*(-drar*antpos[0,i]/2.998e8 - ddecr*antpos[1,i]/2.998e8)

    total_delay_string = '\n'.join(delays.flatten().astype('str'))+'\n'
    with open('delays_burst.dat', 'w', encoding='utf-8') as f:
        f.write(total_delay_string)
    f.close()

    

def set_default_if_unset(start_offset: int, end_offset: int) -> tuple:
    """Set default `start_offset` and `end_offset` if they aren't set."""
    if start_offset < 0:
        start_offset = None
    if end_offset < 0:
        end_offset = None
    return start_offset, end_offset

def get_reftime() -> "Time":
    """Get the reference time used in the real-time correlator."""
    conf = dsc.Conf(use_etcd=True)
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
        '--ntint',
        type=int,
        nargs='?',
        default=32,
        help="number of native time bins to integrate during correlation")
    parser.add_argument(
        '--obsra',
        type=float,
        nargs='?',
        default=20.00,
        help="obsra")
    parser.add_argument(
        '--obsdec',
        type=float,
        nargs='?',
        default=54.57,
        help="obsdec")
    parser.add_argument(
        '--burstra',
        type=float,
        nargs='?',
        default=20.00,
        help="burstra")
    parser.add_argument(
        '--burstdec',
        type=float,
        nargs='?',
        default=54.57,
        help="burstdec")
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
    parser.add_argument(
        '--weights',
        type=str,
        nargs='?',
        default=None,
        help="datestring of weights files to apply pre-correlation")

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    print("Running main program")
    ARGS = parse_commandline_arguments()
    voltages_to_ms(
        ARGS.candname, ntint=ARGS.ntint, start_offset=ARGS.startoffset,
        end_offset=ARGS.stopoffset, dispersion_measure=ARGS.dm, weights=ARGS.weights,obsra=ARGS.obsra,
        obsdec=ARGS.obsdec,burstra=ARGS.burstra,burstdec=ARGS.burstdec)
