"""Simple utilities for T3 imaging.
"""
import subprocess
import json, os, glob
import re
import yaml
from astropy.time import Time
from dsautils import cnf


def load_params(paramfile: str) -> dict:
    """Load parameters for T3 correlation from a yaml file."""
    with open(paramfile) as yamlf:
        T3params = yaml.load(yamlf, Loader=yaml.FullLoader)['T3corr']
    conf = cnf.Conf()
    corrconf = conf.get('corr')
    mfsconf = conf.get('fringe')
    T3params['ch0'] = corrconf['ch0']
    T3params['f0_GHz'] = corrconf['f0_GHz']
    T3params['antennas'] = list(corrconf['antenna_order'].values())[:63]
    T3params['outrigger_delays'] = mfsconf['outrigger_delays']
    return T3params


def get_tstart_from_json(headername: str) -> 'astropy.time.Time':
    """Extract the start time from the header file."""
    with open(headername) as jsonf:
        metadata = json.load(jsonf)
    tstart = Time(metadata['mjds'], format='mjd')
    return tstart


def get_DM_from_json(headername: str) -> float:
    """Extract the dispersion measure from the header file."""
    with open(headername) as jsonf:
        metadata = json.load(jsonf)
    return metadata['dm']


def find_beamformer_weights(candtime: 'astropy.time.Time', bfdir: str = '/data/dsa110/T3/calibs/') -> str:
    """Find the beamformer weights that were in use at a time `candtime`.
    
    In `/data/dsa110/T3/calibs/`, the times in the beamformer weight names are the times when they were
    uploaded to the correlator nodes. Therefore, we want the most recent calibration files that were
    created before `candtime`.
    """
    isot_string = r"[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]"
    isot_pattern = re.compile(isot_string)
    avail_calibs = sorted(
        [
            isot_pattern.findall(calib_path)[0] for calib_path
            in glob.iglob(f"{bfdir}/beamformer_weights_{isot_string}.yaml")],
        reverse=True)
    for avail_calib in avail_calibs:
        if avail_calib < isot_pattern.findall(candtime.isot)[0]:
            return avail_calib


def rsync_file(infile: str, outfile: str) -> str:
    """Rsyncs a file from the correlator machines.

    Parameters
    ----------
    infile : str
        The sourcefile string, e.g. 'corr01.sas.pvt:/home/user/data/fl_out.1.5618974'
    outfile : str
        The destination string, e.g. '/home/user/data/fl_out.1.5618974'

    Returns
    -------
    str
        The full path to the rsynced file in its destination.
    """
    command = f"rsync -avvP --inplace {infile} {outfile}"
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True
    )
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)
    if process.returncode != 0:
        return None
    return outfile
