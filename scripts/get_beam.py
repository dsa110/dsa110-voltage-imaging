#!/usr/bin/env python3                                                          
# -*- coding: utf-8 -*-                                                         
# dsahead python 3.7                                                           

import argparse
import dsacalib.constants as ct
import dsautils.cnf as cnf
import numpy as np
from astropy import time
from astropy import units as u
from astropy import wcs
from astropy.coordinates import ITRS, EarthLocation, SkyCoord
from astropy.time import Time
from dsautils import coordinates, dsa_store
from tqdm import tqdm

# TODO should move to calib constants                                           
Ddish = 4.5
dmax = 273.0
freq = 1428.0
lamb = 3e2 / freq
nbeam = 256
beam_separation = 1  # arcminutes                                               
HA_pointing = 0.0
npx = 1000
theta_min = -(nbeam / 2.0) * beam_separation

def gauss1d(x, mu, sigma):
    """1D Gaussian envelope

    Parameters
    ----------
    x: float
        x-value to evaluate at
    mu: float
        mean
    sigma: float
        sigma

    Returns
    -------
    float
        value of 1d gaussian
    """

    return np.exp(-((x - mu) ** 2) / (2.0 * sigma**2))


def gaussian2D(
    coords,  # x and y coordinates for each image.
    amplitude=1,  # Highest intensity in image.
    xo=0,  # x-coordinate of peak centre.
    yo=0,  # y-coordinate of peak centre.
    sigma_x=1,  # Standard deviation in x.
    sigma_y=1,  # Standard deviation in y.
    rho=0,  # Correlation coefficient.
    offset=0,
    rot=0,
):  # rotation in degrees.

    """Returns 2D Gaussian envelope

    Parameters
    ----------
    coords : (x,y) tuple of meshgrid coordinates for each image
    amplitude: float
        Highest intensity in image. (default 1.0)
    x0: int
        x-coordinate of peak centre. (default 0)
    y0: int
        y-coordinate of peak centre.  (default 0)
    sigma_x: float
        Standard deviation in x. (default 1.0)
    sigma_y: float
        Standard deviation in y. (default 1.0)
    rho: float
        Correlation coefficient between x and y, stretches Gaussian (default 0)
    offset: float
        Baseline offset (default 0.0)
    rot: float
        Rotation of Gaussian in degrees (default 0.0)

    Returns
    -------
    Squeezed 2D grid representing Gaussian.

    """

    x, y = coords

    rot = np.deg2rad(rot)

    x_ = np.cos(rot) * x - y * np.sin(rot)
    y_ = np.sin(rot) * x + np.cos(rot) * y

    xo = float(xo)
    yo = float(yo)

    xo_ = np.cos(rot) * xo - yo * np.sin(rot)
    yo_ = np.sin(rot) * xo + np.cos(rot) * yo

    x, y, xo, yo = x_, y_, xo_, yo_

    # Create covariance matrix
    mat_cov = [
        [sigma_x**2, rho * sigma_x * sigma_y],
        [rho * sigma_x * sigma_y, sigma_y**2],
    ]
    mat_cov = np.asarray(mat_cov)
    # Find its inverse
    mat_cov_inv = np.linalg.inv(mat_cov)

    # PB We stack the coordinates along the last axis
    mat_coords = np.stack((x - xo, y - yo), axis=-1)

    G = (
        amplitude
        * np.exp(
            -0.5
            * np.matmul(
                np.matmul(mat_coords[:, :, np.newaxis, :], mat_cov_inv),
                mat_coords[..., np.newaxis],
            )
        )
        + offset
    )
    return G.squeeze()


# this is in MD/alt coordinates
def get_2Dbeam_model(aliased=False, neighbors=False):
    """Returns a 2D beam model image composed of 2D Gaussians

    Parameters
    ----------
    aliased bool
        For each synth beam shape, add an aliased one offset by +-128 synthesized beams
    neighbors bool
        For each synth beam shape, add an aliased one offset by +-1 synthesized beams
    Returns
    -------
    numpy array
        [nbeam, npix, npix] array of beam response values
    """

    # Arcminutes
    sb_width_fwhm = lamb / dmax * 180 / np.pi * 120

    # Arcminutes
    primary_width_fwhm = lamb / Ddish * 180 / np.pi * 60

    beam_formed_min, beam_formed_max = (
        -2.0 * primary_width_fwhm,
        2.0 * primary_width_fwhm,
    )

    theta = np.linspace(beam_formed_min, beam_formed_max, npx)

    coords = np.meshgrid(theta, theta)

    #    G0 = gaussian2D(coords, xo=0, yo=0, sigma_x=sb_width_fwhm/2.355, sigma_y=primary_width_fwhm/2.355)
    beam_env = gauss1d(theta, 0, primary_width_fwhm / 2.355)
    mus = np.arange(256) * beam_separation + theta_min
    beam_val = np.zeros([nbeam, len(theta), len(theta)])
    #    Genv = gaussian2D(coords, xo=0, yo=0, sigma_x=primary_width_fwhm/2.355, sigma_y=primary_width_fwhm/2.355)

    print("Calculating beam model...")

    for ii in tqdm(range(len((mus)))):
        G = gaussian2D(
            coords,
            xo=mus[ii],
            yo=0,
            sigma_x=sb_width_fwhm / 2.355,
                sigma_y=primary_width_fwhm / 2.355 * 10,
        )
        
        beam_val[ii] = G
        if aliased:
            # pick either +-128 (0->128, 128->0, 255->127)
            iia = np.mod(ii - 128, 256)
            Ga = gaussian2D(
                coords,
                xo=mus[iia],
                yo=0,
                sigma_x=10 * sb_width_fwhm / 2.355,
                sigma_y=primary_width_fwhm / 2.355 * 10,
            )
            beam_val[ii] += Ga
        if neighbors:
            beam_val[ii + 1] += G
            beam_val[ii - 1] += G

    return beam_val


def beams_coord(
    ra_deg,
    dec_deg,
    mjd,
    dec=71.4,
    response=0.1,
    beam_model=None,
    aliased=False,
    neighbors=False,
):
    """Provides a list of beams with response above a given value at a coordinate

    Parameters
    ----------
    ra_deg: float (degrees)
        RA coordinate (ICRS) in degrees
    dec_deg: float (degrees)
        DEC coordinate (ICRS) in degrees
    mjd: float
        mjd to calculate beams at
    dec: float (degrees)
        Pointing declination (default is to try to figure it out from etcd/cnf)
    response: float
        Response limit above which beam is returned, as fraction of peak (default 0.1)
    beam_model: output of get_2Dbeam_model (default None)
    aliased: bool
        Add aliased beam to 2D beam model used for filtering triggers
    neighbors: bool
        Add neighborsing beam to 2D beam model used for filtering triggers


    Returns
    -------
    list
        list of beams that contribute at position (can be empty)
    """

    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame="icrs")
    RA_radians = np.deg2rad(c.ra.value)
    Dec = dec

    # Arcminutes
    sb_width_fwhm = lamb / dmax * 180 / np.pi * 120

    # Arcminutes
    primary_width_fwhm = lamb / Ddish * 180 / np.pi * 60

    beam_formed_min, beam_formed_max = (
        -2.0 * primary_width_fwhm,
        2.0 * primary_width_fwhm,
    )
    theta = np.linspace(beam_formed_min, beam_formed_max, npx)
    coords = np.meshgrid(theta, theta)

    dalt_arr, dMD_arr = coords[0], coords[1]  # these are alt, MD
    if beam_model is not None:
        beam_val = beam_model
    else:
        print("No beam model provided.")
        return [], []

    # get results
    t = Time(mjd, format="mjd", scale="utc")
    c_ITRS = c.transform_to(ITRS(obstime=t))
    local_ha = ct.OVRO_LON * u.rad - c_ITRS.spherical.lon
    HA_src = local_ha.deg + 360.0
    RA_pt = (
        t.sidereal_time(
            "apparent", longitude=ct.OVRO_LON * (180.0 / np.pi) * u.deg
        )
    ).deg
    coord_pt = SkyCoord(ra=RA_pt * u.deg, dec=Dec * u.deg, frame="icrs")
    w = coordinates.create_WCS(coord_pt, ((theta[1] - theta[0]) / 60.0)*u.deg, npx)
    xy = w.wcs_world2pix((RA_pt + HA_src), c.dec.deg, 0)

    if xy[1] >= npx:
        return [], []
    if xy[1] < 0:
        return [], []
    if xy[0] >= npx:
        return [], []
    if xy[0] < 0:
        return [], []

    result = beam_val[:, int(xy[1]), int(xy[0])]
    beams_out = []
    result = result.ravel()
    vals = []
    i = 0
    for val in result:
        if val > response:
            beams_out.append(i)
            vals.append(val)
        i += 1
    return beams_out, vals

def run(mjd,ra,dec,ptdec):

    beams,vals = beams_coord(ra,dec,mjd,dec=ptdec,beam_model=get_2Dbeam_model())

    if len(beams)>0:
        print("Found source...")
        for i,bm in enumerate(beams):
            print(f"Beam {beams[i]}, response {vals[i]}")
    else:
        print("Source not in a beam")


def parse_commandline_arguments() -> 'argparse.Namespace':
    """ Parse command line args. """
    parser = argparse.ArgumentParser(
        description="Get beam from MJD, ra, and dec",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--ra',
        type=float,
        nargs='?',
        default=None,
        help="RA (deg)")
    parser.add_argument(
        '--dec',
        type=float,
        nargs='?',
        default=None,
        help="DEC (deg)")
    parser.add_argument(
        '--mjd',
        type=float,
        nargs='?',
        default=None,
        help="MJD")
    parser.add_argument(
        '--ptdec',
        type=float,
        nargs='?',
        default=71.6,
        help="Pointing dec (deg)")
    
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    print("Running main program")
    ARGS = parse_commandline_arguments()
    run(ARGS.mjd,ARGS.ra,ARGS.dec,ARGS.ptdec)


