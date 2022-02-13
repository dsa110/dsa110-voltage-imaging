"""Calibrate and image measurement sets from T3 visibilities.
"""
import yaml
import numpy as np
from pkg_resources import resource_filename
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.coordinates import Angle
import casatools as cc
from casacore.tables import table
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
            (brightest_point[0]-expected_point[0]).to(u.arcmin),
            (brightest_point[1]-expected_point[1]).to(u.arcmin)
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

def calibrate_T3ms_percorrnode(msnames: dict, bfweights: str, bfdir: str='/data/dsa110/T3/calibs/'):
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

def calibrate_T3ms(msname, bfweights, bfdir: str='/data/dsa110/T3/calibs/'):
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
    gains = gains.reshape(gains.shape[0], -1, gains.shape[-1])
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
        data = data.reshape((nbl*nt, nspw, nfreq//nspw, npol)).swapaxes(0, 1).reshape(nbl*nt*nspw, nfreq//nspw, npol)
        flags = flags.reshape((nbl*nt, nspw, nfreq//nspw, npol)).swapaxes(0, 1).reshape(nbl*nt*nspw, nfreq//nspw, npol)
    
    with table('{0}.ms'.format(msname), readonly=False) as tb:
        tb.putcol('CORRECTED_DATA', data)
        tb.putcol('FLAG', flags)
