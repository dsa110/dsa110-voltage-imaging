import json
import os.path
import numpy as np
import subprocess
import pipes

_allowed = ['astrophysical', 'instrumental', 'unsure/noise', 'rfi', 'archive']

def exists_remote(host, path):
    """Test if a file exists at path on a host accessible with SSH."""
    status = subprocess.call(
        ['ssh', host, 'test -f {}'.format(pipes.quote(path))])
    if status == 0:
        return True
    if status == 1:
        return False
    raise Exception('SSH failed')

def check_voltages(candname):

    filename = f'/home/ubuntu/data/T3/{candname}.json'
    assert os.path.exists(filename), f'candidate json file {filename} not found'
    dd = readfile(filename)
    
    corrs = ['corr03','corr04','corr05','corr06','corr07','corr08','corr10','corr11','corr12','corr14','corr15','corr16','corr18','corr19','corr21','corr22']

    # edit corr03_data and corr03_header
    for corr in corrs:

        data_file = '/home/ubuntu/data/'+candname+'_data.out'
        header_file = '/home/ubuntu/data/'+candname+'_header.json'
        if exists_remote(corr+'.sas.pvt',data_file):
            dd[corr+'_data'] = True
            print('Found data:',corr)
        if exists_remote(corr+'.sas.pvt',header_file):
            dd[corr+'_header'] = True
            print('Found header:',corr)

    writefile(dd, filename)

def readfile(filename):
    """ Read candidate json trigger file and return dict
    Also accepts npy file
    TODO: add file lock?
    """

    try:
        with open(filename, 'r') as fp:
            dd = json.load(fp)
        return dd
    except:
        print('File is not json')
    try:
        dd = np.load(filename,allow_pickle=True)
        return dd.tolist()
    except:
        print('File is not .npy')
        return None
        

def writefile(dd, filename):
    """ Write candidate json trigger file with dict
    """

    with open(filename, 'w') as fp:
        json.dump(dd, fp)


def list_cands_labels(filename):
    """ read json file and list all candidates and labels.
    TODO: decide if more than one allowed
    """
    
    dd = readfile(filename)
    candnames = list(dd.keys())
    for candname in candnames:
        labels = [kk for kk in dd[candname].keys() if kk in _allowed]
        if len(labels):
            labelstr = ', '.join(labels)
        else:
            labelstr = 'no labels'
        print(f'{candname}: {labelstr}')


def set_label(candname, label, filename=None):
    """ Read, add label, and write candidate json file.
    Can optionally provide full path to file.
    Default assumes name of <candname>.json in cwd.
    TODO: decide if file can have more than one candname.
    """

    assert label in _allowed, f'label must be in {_allowed}'

    if filename is None:
        filename = f'/home/ubuntu/data/T3/{candname}.json'

    assert os.path.exists(filename), f'candidate json file {filename} not found'

    dd = readfile(filename)

    if label == 'archive':
        dd['save'] = True
    else:
        dd['label'] = label
    writefile(dd, filename)
        
        
def set_notes(candname, notes, filename=None):
    """ Read, add notes, and write candidate json file.
    Can optionally provide full path to file.
    Default assumes name of <candname>.json in cwd.
    TODO: decide if file can have more than one candname.
    """

    if filename is None:
        filename = f'/home/ubuntu/data/T3/{candname}.json'

    assert os.path.exists(filename), f'candidate json file {filename} not found'

    dd = readfile(filename)
    dd['notes'] = notes
    writefile(dd, filename)
