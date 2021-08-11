from time import sleep
import etcd3
from dask.distributed import Client

client = Client('127.0.0.1:8786')
etcd = etcd3.Etcd3Client()  # TODO: replace with DsaStore

def task(a):
    sleep(5)
    return 'hi, '+str(a)

tasks = []
def cb_func(resp):
    global tasks
    for event in resp.events:
        res = client.submit(task, event.value.decode('utf-8'))
        tasks.append(res)

wid = etcd.add_watch_callback('/mon/corr/1/trigger', cb_func)

while True:
    try:
        print(f'{len(tasks)} tasks in queue')
        for future in tasks:
            if future.done():
                tasks.remove(future)
                print(future.result())
        sleep(1)
    except KeyboardInterrupt:
        print('Exiting after last pass through tasks')
        break
