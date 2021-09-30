"""Simple utilities for T3 imaging.
"""
import subprocess
import json, os, glob


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
#        print('Key <save> not in {0}'.format(fl))
        # Skip corr node json files without the save key if OoD archives twice
        if de.get('save', False):
            print('Will save voltages for ', de['trigname'])
            saved_trigname.append(de['trigname'])

            for corr in ['corr03', 'corr04', 'corr05', 'corr06', 'corr07', 'corr08', 'corr10', 'corr11', 'corr12', 'corr14', 'corr15', 'corr16', 'corr18', 'corr19', 'corr21', 'corr22']:
                if de[corr+'_header'] is True:
                    outfile_h = T3root + datestring + '/'+corr+'_'+de['trigname']+'_header.json'
                    if not os.path.exists(outfile_h):
                        print('copying header '+corr+' '+de['trigname'])
                        os.system('scp '+corr+'.sas.pvt:./data/'+de['trigname']+'_header.json '+outfile_h)
            
                if de[corr+'_data'] is True:
                    outfile_d = T3root + datestring + '/'+corr+'_'+de['trigname']+'_data.out'
                    if not os.path.exists(outfile_d):
                        print('copying data '+corr+' '+de['trigname'])
                        os.system('scp '+corr+'.sas.pvt:./data/'+de['trigname']+'_data.out '+outfile_d)

    return saved_trigname        

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
