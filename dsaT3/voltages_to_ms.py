"""
A script to convert voltage files to measurement sets.
"""
import json
import re
import shutil
import os
import glob
import subprocess
from multiprocessing import Process, Queue
import argparse
import yaml
from pkg_resources import resource_filename
from astropy.time import Time
from dsaT3.utils import get_declination_mjd
from dsaT3.T3imaging import generate_T3_uvh5, calibrate_T3ms
from dsacalib.ms_io import uvh5_to_ms

CORRDIR = '/media/ubuntu/data/dsa110/voltage/'
PROCESSQ = Queue()
NPROC = 3
PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
with open(PARAMFILE) as YAMLF:
    T3PARAMS = yaml.load(YAMLF, Loader=yaml.FullLoader)['T3corr']

def _process_handler(
        deltat_ms,
        deltaf_MHz,
        name,
        declination,
        tstart,
        ntint,
        nfint,
        start_offset,
        end_offset
):
    """Correlates data using T3 cpu correlator.

    Parameters
    ----------
    deltat_ms : float
        The desired integration time in the correlated data, in ms.
    deltaf_MHz : float
        The desired integration frequency in the correlated data, in MHz.
    """
    while not PROCESSQ.empty():
        srcfile, vfile = PROCESSQ.get()
        shutil.copy(
            srcfile,
            vfile
        )
        command = (
            '/home/ubuntu/proj/dsa110-shell/dsa110-bbproc/dsacorr '
            '-d {0} -o {0}.corr -t {1} -f {2} -a 30'.format(
                vfile,
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
        os.remove(vfile)
        corr_files = dict({})
        corr = re.findall('corr\d\d', vfile)[0]
        corr_files[corr] = '{0}.corr'.format(vfile)
        uvh5name = generate_T3_uvh5(
            name,
            declination,
            tstart,
            ntint=ntint,
            nfint=nfint,
            filelist=corr_files,
            start_offset=start_offset,
            end_offset=end_offset
        )
        print(uvh5name)
        os.remove(corr_files[corr])

def __main__(name, filelist, ntint, nfint, start_offset, end_offset):
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
    # Get metadata
    with open('{0}.json'.format(filelist[0])) as jsonf:
        metadata = json.load(jsonf)#[name]
        key = list(metadata.keys())[0]
        metadata = metadata[key]
    tstart = Time(metadata['mjds'], format='mjd')
    declination = get_declination_mjd(tstart)
    deltat_ms = ntint*T3PARAMS['deltat_s']*1e3
    deltaf_MHz = T3PARAMS['deltaf_MHz']

    # Copy files
    # Do 3 at a time
    for file in filelist:
        corr = re.findall('corr\d\d', file)[0]
        fname = file.split('/')[-1]
        outfile = '{0}/{1}_{2}'.format(CORRDIR, corr, fname)
        PROCESSQ.put([file, outfile])
    processes = []
    for i in range(NPROC):
        processes += [Process(
            target=_process_handler,
            args=(
                deltat_ms,
                deltaf_MHz,
                name,
                declination,
                tstart,
                ntint,
                nfint,
                start_offset,
                end_offset
            ),
            daemon=True
        )]
        processes[i].start()
    for thread in processes:
        thread.join()
    hdf5files = sorted(glob.glob('{0}/{1}_corr??.hdf5'.format(
        T3PARAMS['msdir'], name
    )))
    uvh5_to_ms(
        hdf5files,
        '{0}/{1}'.format(T3PARAMS['msdir'], name)
    )

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
        default=8,
        help='number of native time bins to integrate during correlation'
    )
    parser.add_argument(
        '--nfint',
        type=int,
        nargs='?',
        default=8,
        help='number of native freq bins to integrate during correlation'
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
        default=2484,
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
             start_offset=args.startoffset, end_offset=args.stopoffset)
