"""Simple utilities for T3 imaging.
"""
import subprocess



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
