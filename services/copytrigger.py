"""COPPYTRIGGER.PY
Dana Simard, 2021-09-17, dana.simard@astro.caltech.edu

A script to copy over voltage files from corrnodes to the T3 server.
"""

import os
import sys
import time
import datetime
from typing import Callable
from multiprocessing import Process, Queue
from astropy.time import Time
from dsautils import dsa_store
from dsautils import cnf
import dsautils.dsa_syslog as dsl
from dsaT3.utils import rsync_file

CONF = cnf.Conf()
CORR_CONF = CONF.get('corr')
CORRLIST = list(CORR_CONF['ch0'].keys())
NCORR = len(CORRLIST)

# Logger
LOGGER = dsl.DsaSyslogger()
LOGGER.subsystem("software")
LOGGER.app("T3")

ETCD = dsa_store.DsaStore()
DATESTRING = ETCD.get_dict('/cnf/datestring')
RSYNC_Q = Queue()
GATHER_Q = Queue()
ASSESS_Q = Queue()
MAX_ASSESS = 4
MAX_WAIT = 5*60
TSLEEP = 10

T3ROOT = '/media/ubuntu/data/dsa110/T3/'

def change_datestring(newstring: str):
    """Updates the current datestring.

    Parameters
    ----------
    newstring
        The new datestring.
    """
    global DATESTRING
    DATESTRING = newstring
    os.system(f'mkdir -p {T3ROOT}{DATESTRING}')

def copy_voltage_file(payload: dict, queue: Queue = RSYNC_Q):
    """Adds the specified payload to the copy queue.

    Parameters
    ----------
    payload
        The payload to add to the queue.
    queue
        The queue to add the paylod to.
    """
    queue.put(payload)

def task_handler(task_fn: Callable, inqueue: Queue, outqueue: Queue = None, tsleep: int = TSLEEP):
    """Handles in and out queues of preprocessing tasks.

    Parameters
    ----------
    task_fn
        The function to execute, with a single argument.
    inqueue
        The queue containing the arguments to `task_fn`.
    outqueue
        The queue to write the otuput of `task_fn` to.
    """
    while True:
        if not inqueue.empty():
            fname = inqueue.get()
            fname = task_fn(fname)
            if outqueue is not None:
                outqueue.put(fname)
        else:
            time.sleep(tsleep)

def rsync_handler(inqueue: Queue = RSYNC_Q, outqueue: Queue = GATHER_Q):
    """Monitors the rsync queue and rsyncs files.

    Parameters
    ----------
    inqueue
        The queue to monitor
    outqueue
        The queue to update after a successful rsync.
    """
    while True:
        if not inqueue.empty():
            payload = inqueue.get()
            corrname = payload['corrname']
            trigger = payload['trigger']
            candname = next(iter(trigger))
            outheader = rsync_file(
                f'{corrname}:/home/ubuntu/data/{candname}_header.json',
                f'{T3ROOT}{DATESTRING}/{corrname}_{candname}_header.json'
            )
            outdata = rsync_file(
                f'{corrname}:/home/ubuntu/data/{candname}_data.out',
                f'{T3ROOT}{DATESTRING}/{corrname}_{candname}_data.out'
            )
            payload['header_file'] = outheader
            payload['data_file'] = outdata
            outqueue.put(
                payload
            )

def gather_worker(inqueue: Queue, outqueue: Queue, corrlist: list = None):
    """Gathers the payloads from different corrnodes for a single trigger.

    Parameters
    ----------
    inqueue
        The queue to which payloads from the trigger are placed.
    outqueue
        The queue in which to place the gathered payloads.
    corrlist
        The list of corrnode hostnames for which to collect payloads.
    """
    if corrlist is None:
        corrlist = CORRLIST
    nfiles = 0
    end = time.time() + MAX_WAIT
    trigger = None
    candname = None
    ncorr = len(corrlist)
    while nfiles < ncorr and time.time() < end:
        if not inqueue.empty():
            payload = inqueue.get()
            if trigger is None:
                trigger = payload['trigger']
                candname = next(iter(trigger))
                for corrid in corrlist:
                    trigger[candname][f'{corrid}_data'] = None
                    trigger[candname][f'{corrid}_header'] = None
            corrid = payload['corrname']
            trigger[candname][f'{corrid}_data'] = payload['data_file']
            trigger[candname][f'{corrid}_header'] = payload['header_file']
            nfiles += 1
        time.sleep(1)
    outqueue.put(trigger)

def gather_files(inqueue: Queue, outqueue: Queue, ncorr: int = NCORR, max_assess: int = MAX_ASSESS, tsleep: int = TSLEEP):
    """Handles payloads to be gathered and passes them to gathering queues.

    Parameters
    ----------
    inqueue
        The queue to monitor for payloads to be gathered.
    outqueue
        The queue in which to place the gathered payloads.
    ncorr
        The number of corr nodes, and the max number of payloads to be
        gathered together.
    max_assess
        The maximum number of distinct candidates that can be gathered
        at the same time.
    tsleep
        The time to wait if the inqueue is empty.
    """
    gather_queues = [Queue(ncorr) for idx in range(max_assess)]
    gather_names = [None]*max_assess
    gather_processes = [None]*max_assess
    nfiles_assessed = 0
    while True:
        if not inqueue.empty():
            payload = inqueue.get()
            trigger = payload['trigger']
            candname = next(iter(trigger))
            if not candname in gather_names:
                gather_names[nfiles_assessed%max_assess] = candname
                gather_processes[nfiles_assessed%max_assess] = \
                    Process(
                        target=gather_worker,
                        args=(
                            gather_queues[nfiles_assessed%max_assess],
                            outqueue
                        ),
                        daemon=True
                    )
                gather_processes[nfiles_assessed%max_assess].start()
                nfiles_assessed += 1
            gather_queues[
                candname
            ].put([payload])
        else:
            time.sleep(tsleep)

# When updating key for the filterbank writer, use /mon/corr/1/voltagecopy

if __name__ == "__main__":
    processes = {
        'rsync': {
            'nthreads': 1,
            'task_fn': rsync_file,
            'queue': RSYNC_Q,
            'outqueue': GATHER_Q,
            'processes': []
        },
    }
    # Start etcd watch
    os.system(f'mkdir -p {T3ROOT}{DATESTRING}')
    ETCD.add_watch('/mon/corr/1/voltage', copy_voltage_file)
    ETCD.add_watch('/cnf/datestring', change_datestring)
    # Start all threads
    for name in processes.keys():
        for i in range(processes[name]['nthreads']):
            processes[name]['processes'] += [Process(
                target=task_handler,
                args=(
                    processes[name]['task_fn'],
                    processes[name]['queue'],
                    processes[name]['outqueue'],
                ),
                daemon=True
            )]
        for pinst in processes[name]['processes']:
            pinst.start()

    try:
        processes['gather'] = {
            'nthreads': 1,
            'task_fn': gather_files,
            'queue': GATHER_Q,
            'outqueue': ASSESS_Q,
            'processes': []
        }
        processes['gather']['processes'] += [Process(
            target=gather_files,
            args=(
                GATHER_Q,
                ASSESS_Q
                )
        )]
        processes['gather']['processes'][0].start()

        while True:
            for name in processes.keys():
                ETCD.put_dict(
                    '/mon/service/voltagecopy',
                    {
                        "cadence": 60,
                        "time": Time(datetime.datetime.utcnow()).mjd
                        }
                    )
            while not ASSESS_Q.empty():
                final_payload = ASSESS_Q.get()
                ETCD.put_dict(
                    '/mon/corr/1/voltagecopy',
                    final_payload
                )
            time.sleep(60)
    except (KeyboardInterrupt, SystemExit):
        processes['gather']['processes'][0].terminate()
        processes['gather']['processes'][0].join()
        sys.exit()
