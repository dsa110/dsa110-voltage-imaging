from time import sleep
from dask.distributed import Client
from dsautils import dsa_store

#import etcd3
#etcd = etcd3.Etcd3Client()  # TODO: replace with DsaStore

client = Client('127.0.0.1:8786')
de = dsa_store.DsaStore()

def task(a):
    sleep(5)
    return 'hi, '+str(a)

tasks = []
def cb_func(dd):
    global tasks
#    for event in resp.events:
#        res = client.submit(task, event.value.decode('utf-8'))
    res = client.submit(task, dd)
    tasks.append(res)

wid = de.add_watch('/mon/corr/1/trigger', cb_func)

while True:
    try:
        print(f'{len(tasks)} tasks in queue')
        for future in tasks:
            if future.done():
                tasks.remove(future)
                print(future.result())
        sleep(1)
    except KeyboardInterrupt:
        print(f'Cancelling {len(tasks)} tasks and exiting')
        for future in tasks:
            future.cancel()
            tasks.remove(future)
        break
