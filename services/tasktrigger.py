from time import sleep
from dask.distributed import Client
from dsautils import dsa_store
from dsaT3 import T3_manager
import glob, os, json

client = Client('10.42.0.232:8786')
de = dsa_store.DsaStore()

def task(a):

    T3dict = T3_manager.run(a)
    return T3dict

def task_nowait(a):

    T3dict = T3_manager.run_nowait(a)
    return T3dict


tasks = []
def cb_func(dd):
    global tasks    
    res = client.submit(task, dd)
    tasks.append(res)

# set watch
wid = de.add_watch('/mon/corr/1/voltage', cb_func)

# clean up existing triggers
datestring = de.get_dict('/cnf/datestring')
trig_jsons = glob.glob('/data/dsa110/T2/'+datestring+'/cluster_output*.json')
for fl in trig_jsons:
    f = open(fl)
    d = json.load(f)
    trigname = list(d.keys())[0]
    if not os.path.exists('/home/ubuntu/data/T3/'+trigname+'.png'):
        res = client.submit(task_nowait, d)
        tasks.append(res)
    

while True:
    try:
        print(f'{len(tasks)} tasks in queue')
        for future in tasks:
            print(future)
            if future.done():
                print(future.result())
                tasks.remove(future)                
        sleep(5)
    except KeyboardInterrupt:
        print(f'Cancelling {len(tasks)} tasks and exiting')
        for future in tasks:
            future.cancel()
            tasks.remove(future)
        break
