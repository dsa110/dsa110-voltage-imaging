"""
A script to convert voltage files to measurement sets.
"""
import json
import re
import os
import glob
import subprocess
from multiprocessing import Process, Manager
import multiprocessing
import queue
import argparse
import time
import yaml
from pkg_resources import resource_filename
from astropy.time import Time
import astropy.units as u
from dsaT3.utils import get_declination_mjd, rsync_file
from dsaT3.T3imaging import generate_T3_uvh5
from dsacalib.ms_io import uvh5_to_ms

CORRDIR = '/media/ubuntu/ssd/data/B0329/'
NPROC = 8
PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
with open(PARAMFILE) as YAMLF:
    T3PARAMS = yaml.load(YAMLF, Loader=yaml.FullLoader)['T3corr']

#def rsync_handler(
#        rsync_queue,
#        corr_queue,
#        nvoltagefiles,
#        nvoltagefiles_lock,
#        # rsync_done
#):
#    rsync_done = False
#    while not rsync_done:
#        if nvoltagefiles.value >= NPROC:
#            time.sleep(10)
#            continue
#        try:
#            item = rsync_queue.get()
#        except queue.Empty:
#            time.sleep(10)
#        else:
#            if item == 'END':
#                rsync_done = True
#                continue
#            with nvoltagefiles_lock:
#                nvoltagefiles.value += 1
#            srcfile, vfile = item
#            #if not os.path.exists('{0}.corr'.format(vfile)):
#            #    rsync_file(
#            #        srcfile,
#            #        vfile
#            #    )
#            corr_queue.put(vfile)

def corr_handler(
        deltat_ms,
        deltaf_MHz,
        corr_queue,
        uvh5_queue,
        # corr_done,
        #nvoltagefiles,
        #nvoltagefiles_lock,
        ncorrfiles,
        ncorrfiles_lock
):
    """Correlates data using T3 cpu correlator.

    Parameters
    ----------
    deltat_ms : float
        The desired integration time in the correlated data, in ms.
    deltaf_MHz : float
        The desired integration frequency in the correlated data, in MHz.
    """
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
            if not os.path.exists('{0}.corr'.format(vfile)):
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
            #if os.path.exists(vfile):
            #    os.remove(vfile)
            #with nvoltagefiles_lock:
            #    nvoltagefiles.value -= 1
            corr_files = dict({})
            corr = re.findall('corr\d\d', vfile)[0]
            corr_files[corr] = '{0}.corr'.format(vfile)
            uvh5_queue.put(corr_files)

def uvh5_handler(
        name,
        declination,
        tstart,
        ntint,
        nfint,
        start_offset,
        end_offset,
        uvh5_queue,
        #uvh5_done,
        ncorrfiles,
        ncorrfiles_lock
):
    proc = multiprocessing.current_process()
    uvh5_done = False
    while not uvh5_done:
        #print(uvh5_done.is_set())
        try:
            corr_files = uvh5_queue.get()
        except queue.Empty:
            time.sleep(10)
        else:
            if corr_files == 'END':
                print('proc {0} is setting uvh5_done'.format(proc.pid))
                uvh5_done = True
                #uvh5_done.set()
                continue
            uvh5name = generate_T3_uvh5(
                '{0}/{1}'.format(CORRDIR, name),
                declination,
                tstart,
                ntint=ntint,
                nfint=nfint,
                filelist=corr_files,
                start_offset=start_offset,
                end_offset=end_offset
            )
            print(uvh5name)
            #for value in corr_files.values():
            #    os.remove(value)
            with ncorrfiles_lock:
                ncorrfiles.value -= 1
    print('{0} exiting'.format(proc.pid))

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
    if start_offset < 0:
        start_offset = None
    if end_offset < 0:
        end_offset = None
    # Get metadata
    filename = '{0}.json'.format(filelist[0][:-4])
    #outfile = rsync_file(filename, CORRDIR)
    with open(outfile) as jsonf:
        metadata = json.load(jsonf)#[name]
        key = list(metadata.keys())[0]
        metadata = metadata[key]
    tstart = Time(metadata['mjds'], format='mjd')
    #try:
    #    declination = get_declination_mjd(tstart)
    #except ConnectionError:
    declination = 54.58209895*u.deg
    deltat_ms = ntint*T3PARAMS['deltat_s']*1e3
    deltaf_MHz = T3PARAMS['deltaf_MHz']

    manager = Manager()
    #nvoltagefiles = manager.Value('i', 0)
    #nvoltagefiles_lock = manager.Lock()
    ncorrfiles = manager.Value('i', 0)
    ncorrfiles_lock = manager.Lock()
    #rsync_queue = manager.Queue()
    corr_queue = manager.Queue()
    uvh5_queue = manager.Queue()
    #rsync_done = manager.Event()
    #corr_done = manager.Event()
    #uvh5_done = manager.Event()
    # Copy files
    # Do 3 at a time
    for filename in filelist:
        #corr = re.findall('corr\d\d', filename)[0]
        #fname = filename.split('/')[-1]
        #outfile = '{0}/{1}_{2}'.format(CORRDIR, corr, fname)
        #print(filename, outfile)
        #rsync_queue.put([filename, outfile])
        corr_queue.put(filename)
    for i in range(NPROC):
        corr_queue.put('END')
    #rsync_queue.put('END')
    processes = []
    #processes += [Process(
    #    target=rsync_handler,
    #    args=(
    #        rsync_queue,
    #        corr_queue,
    #        nvoltagefiles,
    #        nvoltagefiles_lock,
    #        #rsync_done
    #        ),
    #    daemon=True
    #)]
    for i in range(NPROC):
        processes += [Process(
            target=corr_handler,
            args=(
                deltat_ms,
                deltaf_MHz,
                corr_queue,
                uvh5_queue,
                #corr_done,
                #nvoltagefiles,
                #nvoltagefiles_lock,
                ncorrfiles,
                ncorrfiles_lock
            ),
            daemon=True
        )]
    for i in range(NPROC):
        processes += [Process(
            target=uvh5_handler,
            args=(
                name,
                declination,
                tstart,
                ntint,
                nfint,
                start_offset,
                end_offset,
                uvh5_queue,
                #uvh5_done,
                ncorrfiles,
                ncorrfiles_lock
            ),
            daemon=True
        )]
    for proc in processes:
        proc.start()
    for proc in processes[:NPROC]:
        proc.join()
    print('All corr processes done')
    # We get here
    for i in range(NPROC):
        uvh5_queue.put('END')
    for proc in processes[NPROC:]:
        proc.join()
        print('A uvh5 process joined.')
    print('All uvh5 processes done')
    hdf5files = sorted(glob.glob('{0}/{1}_corr??.hdf5'.format(
        CORRDIR, name
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
