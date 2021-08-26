import numpy as np
from dsautils import dsa_store
from dsaT3 import filplot_funcs as filf
ds = dsa_store.DsaStore()
import time, os
import json

TIMEOUT_FIL = 60
TIMEOUT_CORR = 21600
FILPATH = '/data/dsa110/T1/'
OUTPUT_PATH = '/home/ubuntu/data/T3/'
FIL_CORRS = ['corr01','corr02','corr09','corr13']

# fills output_dict with empty entries
def fill_empty_dict(od):

    od['filfile'] = None
    od['candplot'] = None
    od['save'] = False
    for corr in ['corr03','corr04','corr05','corr06','corr07','corr08','corr10','corr11','corr12','corr14','corr15','corr16','corr18','corr19','corr21','corr22']:
        od[corr+'_data'] = None
        od[corr+'_header'] = None


# searches for local file
def search_for_local_file(fl):

    if os.path.exists(fl):
        return fl
    return None
        

# waits for local file to be written
def wait_for_local_file(fl,timt):

    time_counter = 0
    while not os.path.exists(fl):
        time.sleep(1)
        time_counter += 1
        if time_counter > timt:
            return None

    # wait in case file hasn't been written
    time.sleep(10)

    return fl
    
    
# a is T3 trigger dict
def run(a):

    # set up output dict and datestring
    datestring = ds.get_dict('/cnf/datestring')
    output_dict = a[list(a.keys())[0]]
    output_dict['trigname'] = list(a.keys())[0]
    output_dict['datestring'] = datestring
    fill_empty_dict(output_dict)

    # wait for specific filterbank file to be written
    ibeam = output_dict['ibeam'] + 1
    corrXX = FIL_CORRS[int( (ibeam-1) / 64)]
    filfile = '/data/dsa110/T1/' + corrXX + '/' + datestring + '/fil_' + output_dict['trigname'] + '/' + output_dict['trigname'] +	'_' + str(ibeam) + '.fil'
    print(filfile)
    found_filfile = wait_for_local_file(filfile,TIMEOUT_FIL)

    if found_filfile is None:

        with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'
            json.dump(output_dict, f, ensure_ascii=False, indent=4)
        
        return output_dict
    
    # launch candplotter
    try:
        output_dict['candplot'] = filf.filplot_entry(datestring,a)
    except:
        print('Could not make filplot '+output_dict['trigname'])
        with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'
            json.dump(output_dict, f, ensure_ascii=False, indent=4)

        return output_dict

    # wait for voltage files to be written
    

    # write output_dict to disk
    with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'                  
        json.dump(output_dict, f, ensure_ascii=False, indent=4)

    return output_dict

# a is T3 trigger dict
def run_nowait(a):

    # set up output dict and datestring
    datestring = ds.get_dict('/cnf/datestring')
    output_dict = a[list(a.keys())[0]]
    output_dict['trigname'] = list(a.keys())[0]
    output_dict['datestring'] = datestring
    fill_empty_dict(output_dict)

    # wait for specific filterbank file to be written
    ibeam = output_dict['ibeam'] + 1
    corrXX = FIL_CORRS[int( (ibeam-1) / 64)]
    filfile = '/data/dsa110/T1/' + corrXX + '/' + datestring + '/fil_' + output_dict['trigname'] + '/' + output_dict['trigname'] +	'_' + str(ibeam) + '.fil'
    print(filfile)
    found_filfile = search_for_local_file(filfile)

    if found_filfile is None:

        with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'
            json.dump(output_dict, f, ensure_ascii=False, indent=4)
        
        return output_dict
    
    # launch candplotter
    try:
        output_dict['candplot'] = filf.filplot_entry(datestring,a)
    except:
        print('Could not make filplot '+output_dict['trigname'])
        with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'
            json.dump(output_dict, f, ensure_ascii=False, indent=4)

        return output_dict

    # wait for voltage files to be written
    

    # write output_dict to disk
    with open(OUTPUT_PATH + output_dict['trigname'] + '.json', 'w') as f: #encoding='utf-8'                  
        json.dump(output_dict, f, ensure_ascii=False, indent=4)

    return output_dict

