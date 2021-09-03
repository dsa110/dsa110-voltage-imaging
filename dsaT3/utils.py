"""Simple utilities for T3 imaging.
"""
import subprocess
import numpy as np
from influxdb import DataFrameClient
import astropy.units as u
from astropy.coordinates import Angle
import dsacalib.constants as ct
from dsacalib.utils import direction
import numpy as np
import datetime
from psrqpy import QueryATNF
from astropy.coordinates import SkyCoord, ITRS, EarthLocation
from dsautils import dsa_store
from progress.bar import Bar
import dsacalib.constants as ct
import dsautils.cnf as cnf
ds = dsa_store.DsaStore()

MY_CNF = cnf.Conf()
CORR_CNF = MY_CNF.get('corr')
query = QueryATNF(params=['DM', 'RAJ', 'DECJ', 'S1400', 'PSRJ', 'PSRB'])
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


def get_pointing_mjd(mjd):
    tmjd = Time(mjd, scale='utc', format='mjd')
    ra_mjd =  (tmjd.sidereal_time('apparent', longitude=ct.OVRO_LON*(180./np.pi)*u.deg)).deg
    dec_mjd = get_declination_mjd(tmjd)

    return ra_mjd*u.deg, dec_mjd


def get_pointing_declination(tol=0.25):
    """Gets the pointing declination from the commanded antenna elevations.
    Parameters
    ----------
    tol : float
        The tolerance for discrepancies in the antenna pointing and commanded
        elevations, in degrees.
    Returns
    -------
    astropy quantity
        The pointing declination, in degrees or equivalent.
    """
    commanded_els = np.zeros(len(CORR_CNF['antenna_order']))
    for idx, ant in CORR_CNF['antenna_order'].items():
        try:
            antmc = ds.get_dict('/mon/ant/{0}'.format(ant))
            a1 = np.abs(antmc['ant_el'] - antmc['ant_cmd_el'])
        except:
            a1 = 2.*tol
        if a1 < tol:
            commanded_els[idx] = antmc['ant_cmd_el']
        else:
            commanded_els[idx] = np.nan

    pt_el = np.nanmedian(commanded_els)
    if pt_el is not np.nan:
        pt_dec = ct.OVRO_LAT*u.rad + pt_el*u.deg - 90*u.deg
    else:
        pt_el = CORR_CNF['pt_dec']
    return pt_dec


def get_pointing_now():
    tnow = Time(datetime.datetime.now(), scale='utc')
    ra_now = (tnow.sidereal_time('apparent', longitude=ct.OVRO_LON*(180./np.pi)*u.deg)).deg
    dec_now = get_pointing_declination()
    dec_now = np.rad2deg(dec_now.value)
    
    return ra_now*u.deg, dec_now*u.deg, tnow.mjd


def match_pulsar(RA_mjd, Dec_mjd, thresh_deg=3.5):
    RA_psr, Dec_psr, DM = np.array(query['RAJ']), np.array(query['DECJ']), np.array(query['DM'])
#    print(RA_mjd, Dec_mjd)
    c = SkyCoord(ra=RA_mjd, dec=Dec_mjd)
    catalog = SkyCoord(ra=RA_psr, dec=Dec_psr, unit=(u.h, u.deg))
    
    ra,dec = catalog.data.lon.deg, catalog.data.lat.value
    sep_deg = np.sqrt((ra-RA_mjd.value)**2 + (dec - Dec_mjd.value)**2)
    ind_near = np.where(sep_deg<thresh_deg)[0]
    #idx, d2, d3 = c.match_to_catalog_sky(catalog)

    return ind_near


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


def get_beam_ha(ibeam, beam_sep=1*u.arcmin):
    return beam_sep*(127-ibeam)


def get_beam_ra_dec(obstime, ibeam):
    """Get ra and dec of beam.
    
    Parameters
    ----------
    obstime : astropy time object
        observing time
    ibeam : int
        beam id
    
    Returns
    -------
    tuple
        Ra, Dec in radians
    """
    trigger_dir = direction(
        'HADEC',
        get_beam_ha(ibeam).to_value(u.rad),
        get_declination_mjd(obstime).to_value(u.rad),
        obstime=obstime
    )
    return trigger_dir.J2000()