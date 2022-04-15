"""Applies beamformer weights and images a measurement set.

Requires a json file with details of the candidates.  When read in should look like:
    candidates = {
        '220311aabg': {
            'candnames': ['220311aabg_long'],
            'datestring': '2022_3_11_15_57_30',
            'coordinates': {'ra': None, 'dec': None},
            'bfweights': None,
            'badants': None},
        '220319aaeb': {
            'candnames': ['220319aaeb_long', '220319aaeb'],
            'datestring': 'temp_220319aaeb',
            'coordinates': {'ra': None, 'dec': None},
            'bfweights': None,
            'badants': [31, 105]}}
"""

from typing import List
import os
import shutil
import argparse

import matplotlib
matplotlib.use('Agg')
from astropy.coordinates import Angle
import astropy.units as u

from casatasks import tclean
from dsacalib.flagging import flag_antenna, flag_baselines, flag_rfi, flag_zeros, reset_all_flags

from dsavim.imaging import calibrate_T3ms, plot_image
from dsavim.utils import find_beamformer_weights, get_tstart_from_json


def main(
        candparams: str, candnames_to_process: List[str] = None,
        msdir: str = "/media/ubuntu/data/dsa110/imaging/") -> None:
    """Read and lookup data on the candidates, and process each one.

    Parameters
    ----------
    candparams : str
        The full path to a json file containing information on the candidates.
    candnames_to_process: list
        The candnames to calibrate and image.  If not provided, all candidates in the json file
        will be processed.
    msdir : str
        The directory containing the candidate measurement sets.
    """

    with open(candparams, encoding="utf-8") as f:
        candidates = json.load(f)

    for cand, canddata in candidates.items():
        if not canddata['bfweights']:
            canddata['bfweights'] = find_beamformer_weights(
                get_tstart_from_json(
                    f"/media/ubuntu/ssd/T3/{canddata['datestring']}/{cand}.json"))

    for candidate in candidates.values():
        
        ra = candidate['coordinates']['ra']
        dec = candidate['coordinates']['dec']
        bfweights = candidate['bfweights']
        candnames = candidate['candnames']

        for candname in candnames:
            if not candnames_to_process or candname in candnames_to_process:
                calibrate_and_image(f'{msdir}{candname}', bfweights, ra=ra, dec=dec)


def calibrate_and_image(
        msname: str, bfweights: str, bfdir: str = "/data/dsa110/T3/calibs/", ra: "Angle" = None,
        dec: "Angle" = None, suffix: str = None) -> None:
    """Calibrate and image a measurement set.

    Parameters
    ----------
    msname : str
        The full path to the measurement set.
    bfweights : str
        The name of the beamformer weights to apply,
        e.g. 'J141120+521209_2021-07-13T02:40:26'
    bfdir : str
        The directory containing to the bfweights,
        e.g. '/data/dsa110/T3/calibs/'
    ra : astropy Angle
        The right ascension the source should be localized to.
    dec : astropy Angle
        The declination the source should be localized to.
    suffix : str
        A suffix to add to the image names.
    """
    if msname[-3:] == '.ms':
        msname = msname[:-3]
    if ra is not None:
        ra = Angle(ra)
    if dec is not None:
        dec = Angle(dec)
        
    reset_all_flags(msname)
    flag_baselines(msname, uvrange="<100m")
    flag_zeros(msname)

    for ant in [31, 105]:
        flag_antenna(msname, ant)
    
    flag_rfi(msname)

    calibrate_T3ms(msname, bfweights, bfdir)

    if not ra:

        # First image
        imname = f'{msname}_full'
        if suffix is not None:
            imname = f'{imname}_{suffix}'
        for extension in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
            if os.path.exists(f'{imname}.{extension}'):
                shutil.rmtree(f'{imname}.{extension}')
        tclean(
            f'{msname}.ms',
            datacolumn='corrected',
            imsize=8192,
            cell='3.0arcsec',
            imagename=imname,
        )
        locx, locy = plot_image(
            f'{imname}.image',
            verbose=False,
            show=False,
            outname=imname
        )
        print(f'Initial localization: {locx} {locy}')

        # Second image
        imname = '{0}_zoom1'.format(msname)
        if suffix is not None:
            imname = f'{imname}_{suffix}'
        for extension in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
            if os.path.exists(f'{imname}.{extension}'):
                shutil.rmtree(f'{imname}.{extension}')
        tclean(
            f'{msname}.ms',
            datacolumn='corrected',
            imsize=16384,
            cell='0.2arcsec',
            imagename=imname,
            phasecenter='J2000 {0}deg {1}deg'.format(
                locx.to_value(u.deg),
                locy.to_value(u.deg)
            )
        )
        locx, locy = plot_image(
            f'{imname}.image',
            verbose=False,
            show=False,
            outname=imname
        )
        print(f'First zoom in: {locx} {locy}')

        # Third image
        imname = f'{msname}_zoom2'
        if suffix is not None:
            imname = f'{imname}_{suffix}'
        for extension in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
            if os.path.exists(f'{imname}.{extension}'):
                shutil.rmtree(f'{imname}.{extension}')
        tclean(
            f'{msname}.ms',
            datacolumn='corrected',
            imsize=16384,
            cell='0.2arcsec',
            imagename=imname,
            phasecenter='J2000 {0}deg {1}deg'.format(
                locx.to_value(u.deg),
                locy.to_value(u.deg)
            )
        )
        locx, locy = plot_image(
            f'{imname}.image',
            verbose=False,
            show=False,
            outname=imname,
        )
        print(f'Second zoom in: {locx} {locy}')

        imname = f'{msname}_zoom3'
        if suffix is not None:
            imname = f'{imname}_{suffix}'
        for extension in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
            if os.path.exists(f'{imname}.{extension}'):
                shutil.rmtree(f'{imname}.{extension}')
        tclean(
            f'{msname}.ms',
            datacolumn='corrected',
        imsize=4096,
        cell='0.1arcsec',
            imagename=imname,
            phasecenter='J2000 {0}deg {1}deg'.format(
                locx.to_value(u.deg),
                locy.to_value(u.deg)
            )
        )
        locx, locy = plot_image(
            f'{imname}.image',
            verbose=False,
            show=False,
            outname=imname
        )
        
    else:
        imname = f'{msname}_radec'.format(msname)
        if suffix is not None:
            imname = f'{imname}_{suffix}'
        for extension in ['image', 'model', 'pb', 'psf', 'residual', 'sumwt']:
            if os.path.exists(f'{imname}.{extension}'):
                shutil.rmtree(f'{imname}.{extension}')
        tclean(
            f'{msname}.ms',
            datacolumn='corrected',
        imsize=4096,
        cell='0.1arcsec',
            imagename=imname,
            phasecenter='J2000 {0}deg {1}deg'.format(
                ra.to_value(u.deg),
                dec.to_value(u.deg)
            )
        )
        locx, locy = plot_image(
            f'{imname}.image',
            verbose=False,
            show=False,
            outname=imname,
            expected_point=(ra, dec)
        )
        print('From image around source ra, dec:')
        print('Offset by {0} {1}'.format(
            (locx-ra).to(u.arcsecond),
            (locy-dec).to(u.arcsecond)
        ))

def parse_commandline_arguments() -> 'argparse.Namespace':
    """Parse commandline arguments."""
    parser = argparse.ArgumentParser(
        description="Apply beamformer weights and calibrate voltage measurement sets.")
    parser.add_argument(
        'candparams',
        type=str,
        help="Full path to json file with candidate parameters")
    parser.add_argument(
        '--candnames',
        type=str,
        help="Names of candidates to process",
        nargs='*')

    args = parser.parse_args()
    return args

if __name__ == '__main__':

    ARGS = parse_commandline_arguments()
    main(ARGS.candparams, ARGS.candnames)
