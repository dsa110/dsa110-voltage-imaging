"""Simple utilities for T3 imaging.
"""

from influxdb import DataFrameClient
import astropy.units as u
from astropy.time import Time
import dsacalib.constants as ct
import numpy as np
import subprocess

influx = DataFrameClient('influxdbservice.sas.pvt', 8086, 'root', 'root', 'dsa110')

def get_elevation_mjd(tobs):
    time_ms = int(tobs.unix*1000)
    query = ('SELECT ant_num, ant_el, ant_cmd_el, ant_el_err FROM "antmon" WHERE '
             'time >= {0}ms and time < {1}ms'.format(time_ms-500, time_ms+500))
    el_df = influx.query(query)
    el_df = el_df['antmon']
    el = np.median(el_df[np.abs(el_df['ant_el_err']) < 1.]['ant_cmd_el'])*u.deg
    return el

def get_declination(elevation, latitude=ct.OVRO_LAT*u.rad):
    return (elevation+latitude-90*u.deg).to(u.deg)

def get_declination_mjd(tobs, latitude=ct.OVRO_LAT*u.rad):
    elevation = get_elevation_mjd(tobs)
    return get_declination(elevation)

def rsync_file(infile, outdir):
    """Rsyncs a file from the correlator machines to dsastorage.
    """
    command = '. ~/.keychain/lxd110h23-sh ; rsync -avvP --inplace {0} {1}'.format(infile, outdir)
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        shell=True
    )
    proc_stdout = str(process.communicate()[0].strip())
    print(proc_stdout)
    fname = infile.split('/')[-1]
    return '{0}{1}'.format(outdir, fname)
