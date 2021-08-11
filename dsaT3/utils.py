"""Simple utilities for T3 imaging.
"""
import subprocess
import numpy as np
from influxdb import DataFrameClient
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import Angle
import dsacalib.constants as ct
from dsacalib.utils import direction

influx = DataFrameClient('influxdbservice.sas.pvt', 8086, 'root', 'root', 'dsa110')

def get_elevation_mjd(tobs):
    """Gets the pointing elevation at a time in the past.

    Parameters
    ----------
    tobs : astropy.time.Time object
        The observing time.
    
    Returns
    -------
    astropy Quantity
        The pointing elevation in degrees or equivalent.
    """
    time_ms = int(tobs.unix*1000)
    query = ('SELECT ant_num, ant_el, ant_cmd_el, ant_el_err FROM "antmon" WHERE '
             'time >= {0}ms and time < {1}ms'.format(time_ms-500, time_ms+500))
    el_df = influx.query(query)
    el_df = el_df['antmon']
    el = np.median(el_df[np.abs(el_df['ant_el_err']) < 1.]['ant_cmd_el'])*u.deg
    return el

def get_declination(elevation, latitude=ct.OVRO_LAT*u.rad):
    """Calculates the declination from the elevation.
    
    Parameters
    ----------
    elevation : astropy Quantity
        The elevation, in degrees or equivalent.
    latitude : astropy Quantity
        The latitude of the telescope, in degrees or equivalent.

    Returns
    -------
    astropy Quantity
        The declination, in degrees or equivalent.
    """
    return (elevation+latitude-90*u.deg).to(u.deg)

def get_declination_mjd(tobs, latitude=ct.OVRO_LAT*u.rad):
    """Gets the pointing declination at a time in the past.

    Parameters
    ----------
    tobs : astropy.Time object
        The observing time.
    latitude : astropy Quantity
        The telescope latitude.

    Returns
    -------
    astropy Quantity
        The declination, in degrees or equivalent.
    """
    elevation = get_elevation_mjd(tobs)
    return get_declination(elevation)

def rsync_file(infile, outdir):
    """Rsyncs a file from the correlator machines to dsastorage.

    Parameters
    ----------
    infile : str
        The sourcefile string, e.g. 'corr01.sas.pvt:/home/user/data/fl_out.1.5618974'
    outfile : str
        The destination string, e.g. '/home/user/data/'
    
    Returns
    -------
    str
        The full path to the rsynced file in its destination.
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

def get_pointing(obstime):
    """Get the RA and DEC of the array at a given time in the past.

    Parameters
    ----------
    obstime : astropy.time.Time object
        The observation time.

    Returns
    -------
    tuple
        (ra, dec) of the observation in J2000 epoch, as astropy Quantities.
    """
    ra, dec = direction(
        'HADEC',
        0.,
        get_declination_mjd(obstime).to_value(u.rad),
        obstime=obstime.mjd
    ).J2000()
    return Angle(ra*u.rad).to(u.hourangle), Angle(dec*u.rad).to(u.deg)