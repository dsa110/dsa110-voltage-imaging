"""
A script to convert voltage files to measurement sets.
"""
import json
import re
import shutil
import os
import subprocess
from multiprocessing import Process, Queue
import argparse
import yaml
from pkg_resources import resource_filename
from astropy.time import Time
from dsaT3.utils import get_declination_mjd
from dsaT3.T3imaging import generate_T3_ms, calibrate_T3ms

CORRDIR = '/home/ubuntu/data/'
CORRQ = Queue()
COPYQ = Queue()
NCORRPROC = 8
NCOPYPROC = 8
MSDIR = '/media/ubuntu/data/dsa110/imaging/'
PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
with open(PARAMFILE) as YAMLF:
    T3PARAMS = yaml.load(YAMLF, Loader=yaml.FullLoader)['T3corr']

def _copy_handler():
    while not COPYQ.empty():
        args = COPYQ.get()
        shutil.copy(
            args[0],
            args[1]
        )

def _corr_handler(deltat_ms, deltaf_MHz):
    """Correlates data using T3 cpu correlator.

    Parameters
    ----------
    deltat_ms : float
        The desired integration time in the correlated data, in ms.
    deltaf_MHz : float
        The desired integration frequency in the correlated data, in MHz.
    """
    while not CORRQ.empty():
        vf = CORRQ.get()
        command = (
            '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/dsacorr '
            '-d {0} -o {0}.corr -t {1} -f {2} -a 30'.format(
                vf,
                deltat_ms,
                deltaf_MHz
            )
        )
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=True
        )
        proc_stdout = str(process.communicate()[0].strip())
        print(proc_stdout)

def __main__(name, filelist, ntint=96, nfint=1, start_offset=158, end_offset=167, clean=True):
    """
    Correlate voltage files and convert to a measurement set.

    Parameters
    ----------
    name : str
        The unique name of the candidate.
    filelist : list
        The full paths to the voltage files on dsa-storage.
    ntint : int
        The number of time samples to integrate together during correlation.
    nfint : int
        The number of frequency channels to integrate together during
        correlation.
    start_offset : int
        The number of time samples (after correlation) to offset the start of
        the measurement set by.  If not provided, the entire time is converted
        to a measurement set.
    end_offset : int
        The last time sample (after correlation) to write to the measurement
        set. If not provide,d the entire time is converted to a measurement
        set.
    """
    # Get metadata
    with open('{0}.json'.format(filelist[0])) as jsonf:
        metadata = json.load(jsonf)#[name]
        key = list(metadata.keys())[0]
        metadata = metadata[key]
    tstart = Time(metadata['mjds'], format='mjd')
    declination = get_declination_mjd(tstart)
    deltat_ms = ntint*T3PARAMS['deltat_s']*1e3
    deltaf_MHz = nfint*T3PARAMS['deltaf_MHz']

    # Copy files
    voltage_files = []
    for file in filelist:
        corr = re.findall('corr\d\d', file)[0]
        fname = file.split('/')[-1]
        outfile = '{0}/{1}_{2}'.format(CORRDIR, corr, fname)
        COPYQ.put([file, outfile])
        voltage_files += [outfile]
    processes = []
    for i in range(NCOPYPROC):
        processes += [Process(
            target=_copy_handler,
            daemon=True
        )]
        processes[i].start()
    for thread in processes:
        thread.join()

    # Correlate files
    for vf in voltage_files:
        CORRQ.put(vf)
    # Spawn 4 processes to handle the correlation
    processes = []
    for i in range(NCORRPROC):
        processes += [Process(
            target=_corr_handler,
            args=(deltat_ms, deltaf_MHz),
            daemon=True
        )]
        processes[i].start()
    for thread in processes:
        thread.join()
        
    corr_files = dict({})
    for vf in voltage_files:
        corr = re.findall('corr\d\d', vf)[0]
        corr_files[corr] = '{0}.corr'.format(vf)
    if clean:
        for vf in voltage_files:
            os.remove(vf)

    # Generate the measurement set
    msname = generate_T3_ms(
        name,
        declination,
        tstart,
        ntint=ntint,
        nfint=nfint,
        filelist=corr_files,
        start_offset=start_offset,
        end_offset=end_offset
    )
    print('{0} created'.format(msname))
    if clean:
        for cf in corr_files.values:
            os.remove(cf)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Correlate candidate voltage files.'
    )
    parser.add_argument(
        'name',
        type=str,
        help='unique candidate name'
    )
    parser.add_argument(
        '--ntint',
        type=int,
        nargs='?',
        const=96,
        help='number of native time bins to integrate during correlation'
    )
    parser.add_argument(
        '--nfint',
        type=int,
        nargs='?',
        const=1,
        help='number of native freq bins to integrate during correlation'
    )
    parser.add_argument(
        '--startoffset',
        type=int,
        nargs='?',
        const=158,
        help='nbins from beginning of correlated data to start writing to ms'
    )
    parser.add_argument(
        '--stopoffset',
        type=int,
        nargs='?',
        const=167,
        help='number of bins from end of correlation to write to ms'
    )
    parser.add_argument(
        'filelist',
        type=str,
        nargs='+',
        help='candidate voltage files'
    )
    args = parser.parse_args()
    __main__(args.name, args.filelist, ntint=args.ntint, nfint=args.nfint,
             start_offset=args.startoffset, stop_offset=args.stop_offset)
