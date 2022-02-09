"""Simple utilities for T3 imaging.
"""
import subprocess
import json, os, glob
import pipes
import yaml
from astropy.time import Time
from dsautils import cnf

def delete_remote(host, path):
    """Delete remote file"""
    status = subprocess.call(['ssh', host, 'rm -rf {}'.format(pipes.quote(path))])
    return status
    raise Exception('SSH failed')

    
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

def load_params(paramfile):
    """Load parameters for T3 correlation from a yaml file."""
    with open(paramfile) as yamlf:
        T3params = yaml.load(yamlf, Loader=yaml.FullLoader)['T3corr']
    conf = cnf.Conf()
    corrconf = conf.get('corr')
    mfsconf = conf.get('fringe')
    T3params['ch0'] = corrconf['ch0']
    T3params['f0_GHz'] = corrconf['f0_GHz']
    T3params['antennas'] = list(corrconf['antenna_order'].values())[:63]
    T3params['outrigger_delays'] = mfsconf['outrigger_delays']
    return T3params

def get_tstart_from_json(headername: str) -> "astropy.time.Time":
    """Extract the start time from the header file."""
    with open(headername) as jsonf:
        metadata = json.load(jsonf)
    tstart = Time(metadata['mjds'], format='mjd')
    return tstart

def find_beamformer_weights(candtime: "astropy.time.Time", bfdir: str = '/data/dsa110/T3/calibs/') -> str:
    """Find the beamformer weights that were in use at a time `candtime`.
    
    In `/data/dsa110/T3/calibs/`, the times in the beamformer weight names are the times when they were
    uploaded to the correlator nodes. Therefore, we want the most recent calibration files that were
    created before `candtime`.
    """
    isot_string = '[0-9][0-9][0-9][0-9]-[0-9][0-9]-[0-9][0-9]T[0-9][0-9]:[0-9][0-9]:[0-9][0-9]'
    isot_pattern = re.compile(isot_string)
    avail_calibs = sorted([isot_pattern.findall(calib_path)[0]
                           for calib_path in glob.iglob(f'{bfdir}/beamformer_weights_{isot_string}.yaml')],
                          reverse=True)
    for avail_calib in avail_calibs:
        if avail_calib < isot_pattern.match(candtime.isot):
            return avail_calib

def rsync_file(infile, outfile):
    """Rsyncs a file from the correlator machines.

    Parameters
    ----------
    infile : str
        The sourcefile string, e.g. 'corr01.sas.pvt:/home/user/data/fl_out.1.5618974'
    outfile : str
        The destination string, e.g. '/home/user/data/fl_out.1.5618974'

    Returns
    -------
    str
        The full path to the rsynced file in its destination.
    """
    command = 'rsync -avvP --inplace {0} {1}'.format(infile, outfile)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True
    )
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)
    if process.returncode != 0:
        return None
    return outfile


def archive(datestring, T3root='/media/ubuntu/data/dsa110/T3/'):
    """ Archive data from corr nodes to local (h23).
    Function form of original services/archive.py script.
    """

    # copy from T3 directory 
    os.system('mkdir -p '+T3root+datestring)
    os.system('rsync -av /home/ubuntu/data/T3/ '+T3root+datestring)

    fls = glob.glob(T3root+datestring+'/*.json')
    saved_trigname = []

    for fl in fls:

        f = open(fl)
        de = json.load(f)

        # if save is true, scp voltage files
        if de['save'] is True:

            saved_trigname.append(de['trigname'])

            for corr in ['corr03', 'corr04', 'corr05', 'corr06', 'corr07', 'corr08', 'corr10', 'corr11', 'corr12', 'corr14', 'corr15', 'corr16', 'corr18', 'corr19', 'corr21', 'corr22']:
                
                outfile_h = T3root + datestring + '/'+corr+'_'+de['trigname']+'_header.json'
                infile_h = corr+'.sas.pvt:/home/ubuntu/data/'+de['trigname']+'_header.json'
                outfile_d = T3root + datestring + '/'+corr+'_'+de['trigname']+'_data.out'
                infile_d = corr+'.sas.pvt:/home/ubuntu/data/'+de['trigname']+'_data.out'

                print('Searching for remote file '+'/home/ubuntu/data/'+de['trigname']+'_data.out')
                if exists_remote(corr+'.sas.pvt','/home/ubuntu/data/'+de['trigname']+'_data.out'):
                    if not os.path.exists(outfile_h):
                        rsync_file(infile_h,outfile_h)
                    if not os.path.exists(outfile_d):
                        print('Copying...')
                        rsync_file(infile_d,outfile_d)                                    
            
    return saved_trigname        

