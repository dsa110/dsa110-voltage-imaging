import dsautils.dsa_store as ds
import sys
import numpy as np
from time import sleep
from dsautils import dsa_functions36
from pkg_resources import resource_filename
import os
import subprocess

# defaults
datestring = 'dummy'
docopy = False

def datestring_func():
    def a(event):
        global datestring
        datestring=event
    return a
        
def docopy_func():
    def a(event):
        global docopy
        if event=='True':
            docopy=True
        if event=='False':
            docopy=False
    return a
    
# add callbacks from etcd
my_ds = ds.DsaStore()
docopy = my_ds.get_dict('/cmd/corr/docopy') == 'True'
datestring = my_ds.get_dict('/cnf/datestring')
my_ds.add_watch('/cnf/datestring', datestring_func())
my_ds.add_watch('/cmd/corr/docopy', docopy_func())
scriptname = resource_filename('dsaT3', 'services/send_cands.bash')

print(datestring)
print(docopy)

while True:
    if docopy:
        cmd = [scriptname, datestring]
        output = subprocess.check_output(
            cmd,
            stderr=subprocess.STDOUT
        )
        print(output)

    key = '/mon/service/send_cands'
    value = {'cadence': 20, 'time': dsa_functions36.current_mjd()}
    try:
        my_ds.put_dict(key, value)
    except:
        print('COULD NOT CONNECT TO ETCD')

    sleep(20)
