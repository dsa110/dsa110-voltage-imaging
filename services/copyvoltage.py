from time import sleep
from dask.distributed import Client
from dsautils import dsa_store
from dsaT3 import T3_manager
import glob, os, json
from dsautils import dsa_functions36
import numpy as np

client = Client('10.42.0.232:8786')
de = dsa_store.DsaStore()

def task(a):

    T3dict = T3_manager.copy(a)
    return T3dict

tasks = []
def cb_func(dd):
    global tasks    
    res = client.submit(task, dd)
    tasks.append(res)
    
# set watches on each corr node to copy data
for i in np.arange(1,17):
    wid = de.add_watch('/mon/corr/'+str(i)+'/voltage', cb_func)    
    
while True:
    try:
        print(f'{len(tasks)} tasks in queue')
        for future in tasks:
            print(future)
            if future.done():
                print(future.result())
                tasks.remove(future)

        de.put_dict('/mon/service/T3copy',{'cadence': 5, 'time': dsa_functions36.current_mjd()})
        sleep(5)
    except KeyboardInterrupt:
        print(f'Cancelling {len(tasks)} tasks and exiting')
        for future in tasks:
            future.cancel()
            tasks.remove(future)
        break
