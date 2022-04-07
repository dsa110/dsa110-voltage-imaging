"""Calibrate and image measurement sets from T3 visibilities.
"""
import os
import yaml
import numpy as np
from pkg_resources import resource_filename
import requests
from PIL import Image
from io import BytesIO

import numpy
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import Angle
from astropy.wcs import WCS
from astropy.table import Table
from astropy.io import fits
from astropy.visualization import PercentileInterval, AsinhStretch
import casatools as cc
from casacore.tables import table
from casatasks import exportfits
from dsacalib.ms_io import extract_vis_from_ms

from dsaT3.utils import load_params

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
T3PARAMS = load_params(PARAMFILE)
CORR_ORDER = [int(k[4:]) for k in list(T3PARAMS['ch0'].keys())]

def plot_image(imname, verbose=False, outname=None, show=True, expected_point=None):
    """Plots an image from the casa-generated image file.

    Parameters
    ---------
    imname : str
        The name full path of the image file.
    verbose : bool
        If set to True, prints some information about the image.
    outname : str
        If provided, saves the image in <outname>_image.png.
    show : bool
        If False, the image is closed at the end of the function.
    cellsize : str
        The size of each pixel, in a Casa-recognized angle.
    """
    error = 0
    ia = cc.image()
    error += not ia.open(imname)
    dd = ia.summary()
    # dd has shape npixx, npixy, nch, npol
    npixx = dd['shape'][0]
    if verbose:
        print('Image shape: {0}'.format(dd['shape']))
    imvals = ia.getchunk(0, int(npixx))[:, :, 0, 0]
    #imvals = fftshift(imvals)
    error += ia.done()
    max_idxs = np.unravel_index(imvals.argmax(), imvals.shape)
    cellsizex = Angle(dd['incr'][0], dd['axisunits'][0])
    cellsizey = Angle(dd['incr'][1], dd['axisunits'][1])
    ra, dec = (
        Angle('{0}{1}'.format(dd['refval'][0], dd['axisunits'][0])),
        Angle('{0}{1}'.format(dd['refval'][1], dd['axisunits'][1]))
    )
    brightest_point = (
        ra +
        Angle('{0}{1}'.format(
            dd['incr'][0]*(max_idxs[0]-dd['refpix'][0]),
            dd['axisunits'][0]
        ))/np.cos(dec),
        dec +
        Angle('{0}{1}'.format(
            dd['incr'][1]*(max_idxs[1]-dd['refpix'][1]),
            dd['axisunits'][1]
        ))
    )
    if verbose:
        print('Peak SNR at pix ({0},{1}) = {2}'.format(max_idxs[0],
                                                       max_idxs[1],
                                                       imvals.max()/
                                                       imvals.std()))
        print('Value at peak: {0}'.format(imvals.max()))
        print('Value at origin: {0}'.format(imvals[imvals.shape[0]//2,
                                                   imvals.shape[1]//2]))

    _, ax = plt.subplots(1, 1, figsize=(15, 8))
    pim = ax.imshow(
        imvals.transpose(),
        interpolation='none',
        origin='lower',
        extent=[
            (-imvals.shape[0]/2*Angle(cellsizex)).to_value(u.arcmin),
            (imvals.shape[0]/2*Angle(cellsizex)).to_value(u.arcmin),
            (-imvals.shape[1]/2*Angle(cellsizey)).to_value(u.arcmin),
            (imvals.shape[1]/2*Angle(cellsizey)).to_value(u.arcmin)
        ]
    )
    plt.colorbar(pim)
    ax.axvline(0, color='white', alpha=0.5)
    ax.axhline(0, color='white', alpha=0.5)
    ax.set_xlabel('l (arcmin)')
    ax.set_ylabel('m (arcmin)')
    plttitle = '{0} {1} {2}'.format(
        imname,
        brightest_point[0].to_string('hourangle'),
        brightest_point[1].to_string('deg')
    )
    if expected_point is not None:
        plttitle += ', offset by {0:.2f} {1:.2f}'.format(
            (brightest_point[0]-expected_point[0]).to(u.arcsecond),
            (brightest_point[1]-expected_point[1]).to(u.arcsecond)
        )
    plt.title(plttitle)
    if outname is not None:
        plt.savefig('{0}_image.png'.format(outname))
    if not show:
        plt.close()
    if error > 0:
        print('{0} errors occured during imaging'.format(error))
    return brightest_point

def read_bfweights(bfweights, bfdir):
    """Reads the beamforming weights.

    Parameters
    ----------
    bfweights : str
        The label of the file containing the weights. Will open
        <bfdir>/beamformer_weights_<bfweights>.yaml
    bfdir : str
        The directory in which the beamformer weights are stored.

    Returns
    -------
    antenna_order : list
        The order of the antennas in the bfweights array.
    bfweights : ndarray
        The beamformer weights, (antenna, freqeuncy, polarization).
        Frequency is in the same order as in the correlator.
    """
    with open('{0}/beamformer_weights_{1}.yaml'.format(
            bfdir,
            bfweights,
    )) as yamlf:
        bfparams = yaml.load(yamlf, Loader=yaml.FullLoader)
    if 'cal_solutions' in bfparams.keys():
        bfparams = bfparams['cal_solutions']
    antenna_order = bfparams.get('antenna_order', T3PARAMS['antennas'])
    corr_order = bfparams.get('corr_order', CORR_ORDER)
    gains = np.zeros(
        (len(antenna_order), len(corr_order), 48, 2),
        dtype=np.complex
    )
    for corridx, corr in enumerate(corr_order):
        with open(
                '{0}/beamformer_weights_corr{1:02d}_{2}.dat'.format(
                    bfdir,
                    corr,
                    bfweights
                ),
                'rb'
        ) as f:
            data = np.fromfile(f, '<f4')
        temp = data[64:].reshape(64, 48, 2, 2)
        gains[:, corridx, :, :] = temp[..., 0]+1.0j*temp[..., 1]

    return antenna_order, gains

def calibrate_T3ms_percorrnode(msnames: dict, bfweights: str, bfdir: str = '/data/dsa110/T3/calibs/'):
    """Calibrates a measurement set using the beamformer weights.

    Calibrated data is written into the CORRECTED_DATA column.

    Parameters
    ----------
    msname : str
        The name of the measurement set.
    bfweights : str
        The label of the file containing the weights. Will open
        <bfdir>/beamformer_weights_<bfweights>.yaml
    bfdir : str
        The directory in which the beamformer weights are stored.
    """
    antenna_order, gains = read_bfweights(bfweights, bfdir)

    for corrnode, msname in msnames.items():
        corridx = CORR_ORDER.index(int(corrnode.strip('corr')))
        apply_calibration(msname, gains[:, corridx, :, :], antenna_order)

def calibrate_T3ms(msname, bfweights, bfdir: str = '/data/dsa110/T3/calibs/'):
    """Calibrates a measurement set using the beamformer weights.

    Calibrated data is written into the CORRECTED_DATA column.

    Parameters
    ----------
    msname : str
        The name of the measurement set.
    bfweights : str
        The label of the file containing the weights. Will open
        <bfdir>/beamformer_weights_<bfweights>.yaml
    bfdir : str
        The directory in which the beamformer weights are stored.
    """
    antenna_order, gains = read_bfweights(bfweights, bfdir)
    gains = gains.reshape((gains.shape[0], -1, gains.shape[-1]))
    apply_calibration(msname, gains, antenna_order)

def apply_calibration(msname: str, gains: np.ndarray, antenna_order: list):
    """Apply `gains` to the visibilities in `msname`."""
    data, _, fobs, flags, ant1, ant2, _, spw, orig_shape = extract_vis_from_ms(
        msname, data='data')
    assert orig_shape in [['time', 'baseline', 'spw'], ['spw', 'time', 'baseline']]
    nspw = len(spw)
    nfreq = len(fobs)
    npol = data.shape[-1]
    nbl = data.shape[0]
    nt = data.shape[1]

    data = data.reshape(nbl, nt, nfreq, npol)
    flags = flags.reshape(nbl, nt, nfreq, npol)

    if np.mean(np.diff(fobs)) > 0:
        assert np.all(np.diff(fobs) > 0)
        gains = gains[:, ::-1, :]
    else:
        assert np.all(np.diff(fobs) < 0)

    for i in range(nbl):
        a1 = ant1[i]+1
        a2 = ant2[i]+1
        try:
            bl_gains = (
                np.conjugate(gains[antenna_order.index(a2), ...])*
                gains[antenna_order.index(a1), ...])
            bl_gains = np.exp(1.j*np.angle(bl_gains))
            data[i, ...] *= bl_gains[np.newaxis, :, :]
        except ValueError:
            flags[i, ...] = 1
            print('no calibration solutions for baseline {0}-{1}'.format(a1, a2))

    data = data.swapaxes(0, 1).reshape((-1, nfreq, npol))
    flags = flags.swapaxes(0, 1).reshape((-1, nfreq, npol))

    if orig_shape == ['spw', 'time', 'baseline']:
        data = data.reshape((nbl*nt, nspw, nfreq//nspw, npol)
                           ).swapaxes(0, 1).reshape(nbl*nt*nspw, nfreq//nspw, npol)
        flags = flags.reshape((nbl*nt, nspw, nfreq//nspw, npol)
                             ).swapaxes(0, 1).reshape(nbl*nt*nspw, nfreq//nspw, npol)

    with table('{0}.ms'.format(msname), readonly=False) as tb:
        tb.putcol('CORRECTED_DATA', data)
        tb.putcol('FLAG', flags)

def get_ps1_images(ra: float, dec: float, size: int = 240, filters: str = "grizy") -> "Table":
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = ("{service}?ra={ra}&dec={dec}&size={size}&format=fits"
           "&filters={filters}").format(**locals())
    table = Table.read(url, format='ascii')
    return table


def get_ps1_url(
        ra: float, dec: float, size: int = 240, output_size: int = None, filters: str = "grizy",
        format: str = "jpg", color: bool = False) -> str:
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    Returns a string with the URL
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = get_ps1_images(ra,dec,size=size,filters=filters)
    url = ("https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
           "ra={ra}&dec={dec}&size={size}&format={format}").format(**locals())
    if output_size:
        url = url + "&output_size={}".format(output_size)
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            url = url + "&{}={}".format(param,table['filename'][i])
    else:
        urlbase = url + "&red="
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    return url


def get_ps1_colorim(
        ra: float, dec: float, size: int = 240, output_size: int = None, filters: str = "grizy",
        format: str = "jpg") -> "image":
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg", "png"):
        raise ValueError("format must be jpg or png")
    url = get_ps1_url(
        ra, dec, size=size, filters=filters, output_size=output_size, format=format, color=True)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def get_ps1_grayim(
        ra: float, dec: float, size: int = 240, output_size: int = None, filter: str = "g",
        format: str = "jpg") -> "image":
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png")
    Returns the image
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = get_ps1_url(ra,dec,size=size,filters=filter,output_size=output_size,format=format)
    r = requests.get(url[0])
    im = Image.open(BytesIO(r.content))
    return im

def load_radio_image(imname: str) -> tuple:
    """Extract a radio image from a fits image."""
    with fits.open(imname) as fh:
        hdu = fh[0]
        wcs = WCS(hdu.header)
        image = hdu.data

    return image[0, 0, :, :], wcs

def load_radio_casa_image(imname: str) -> "np.ndarray":
    ia = cc.image()
    ia.open(imname)
    dd = ia.summary()
    npixx = dd['shape'][0]
    imvals = ia.getchunk(0, int(npixx))[:, :, 0, 0]
    return imvals

def get_centre_coordinates(imname: str) -> tuple:
    """Extract centre coordinates and size of a casa image reconstruction.

    returns (ra_centre, dec_centre, size), all astropy quantities
    """
    ia = cc.image()
    ia.open(imname)
    dd = ia.summary()
    ia.close()
    npixx = dd['shape'][0]
    raref, decref = (
        Angle('{0}{1}'.format(dd['refval'][0], dd['axisunits'][0])),
        Angle('{0}{1}'.format(dd['refval'][1], dd['axisunits'][1]))
    )
    ralist, declist = (
        raref +
        Angle(f"{dd['incr'][0]}{dd['axisunits'][0]}")*(np.arange(npixx)-dd['refpix'][0])/np.cos(decref),
        decref +
        Angle(f"{dd['incr'][1]}{dd['axisunits'][1]}")*(np.arange(npixx)-dd['refpix'][1]),
    )
    ra_centre, dec_centre = np.mean(ralist[npixx//2-1:npixx//2+1]), np.mean(declist[npixx//2-1:npixx//2+1])
    size = Angle(f"{dd['incr'][1]}{dd['axisunits'][1]}")*npixx

    return ra_centre, dec_centre, size

def get_fits_image(ra_centre: "Quantitiy", dec_centre: "Quantity", size: "Quantity") -> "np.ndarray":
    """Extract the ps1 fits image from the online server."""
    fitsurl = get_ps1_url(
    ra_centre.to_value(u.deg), dec_centre.to_value(u.deg), size=int(round(size.to_value(u.arcsecond)*4)), filters="i", format="fits")

    with fits.open(fitsurl[0]) as fh:
        hdu = fh[0]
        wcs = WCS(hdu.header)
        image = hdu.data

    image[np.isnan(image)] = 0.0
    transform = AsinhStretch() + PercentileInterval(99.5)
    image = transform(image)
    return image, wcs

def plot_contours_on_ps1(imname: str, resize_factor: float = None) -> None:
    """Plot radio contours over a ps1 image of the same region."""
    # TODO: get image coordinates from the fits file so we don't need the image file too
    fits_imname = imname.replace(".image", ".fits")
    if not os.path.exists(fits_imname):
        exportfits(imname, fits_imname)

    radio_image, radiowcs = load_radio_image(fits_imname)
    if np.isnan(np.nanmean(radio_image)):
        radio_image = load_radio_casa_image(imname)
        radiowcs = None
    ra_centre, dec_centre, size = get_centre_coordinates(imname)
    ps1_image, wcs = get_fits_image(ra_centre, dec_centre, size)

    immax = np.nanmax(radio_image)
    contours = [0.01, 0.1, 0.25, 0.5, 0.75, 0.90, 0.99]

    fig, ax = plt.subplots(1, 1, figsize=(12, 6), subplot_kw={"projection": radiowcs, "slices": ('x', 'y', 0, 0)})
    ax.imshow(ps1_image, origin="lower", transform=ax.get_transform(wcs))
    ax.contour(radio_image, colors="white", levels=[immax*contour for contour in contours])
    if resize_factor:
        x1, x2 = ax.get_xlim()
        y1, y2 = ax.get_ylim()
        ax.set_xlim((x1+x2)/2-(x2-x1)*resize_factor, (x1+x2)/2+(x2-x1)*resize_factor)
        ax.set_ylim((y1+y2)/2-(y2-y1)*resize_factor, (y1+y2)/2+(y2-y1)*resize_factor)
