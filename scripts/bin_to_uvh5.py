"""Convert a binary file output by the cpu correlator to uvh5."""

import astropy.units as u
from dsautils.coordinates import get_declination, get_elevation
from dsaT3.generate_uvh5 import generate_uvh5
from dsaT3.utils import get_tstart_from_json

def bin_to_uvh5() -> None:
    """Convert a binary file output by the cpu correlator to uvh5."""
    params = get_parameters()

    _ = generate_uvh5(
        f'{params["output_dir"]}/{params["candname"]}',
        params['declination']*u.deg,
        params['tstart'],
        ntint=params['ntint'],
        nfint=params['nfint'],
        filelist=params['corr_files'],
        start_offset=params['start_offset'],
        end_offset=params['end_offset'],
        npol_out=params['npol_out']
    )

def get_parameters() -> tuple:
    """Define parameters for conversion."""
    headername = '/media/ubuntu/data/dsa110/T3/2022_1_21_5_1_21/220121aaat.json'
    tstart = get_tstart_from_json(headername)

    declination = get_declination(get_elevation(tstart)).to_value(u.deg)

    params = {
        'tstart': tstart,
        'declination': declination,
        'corr_files': {'corr22': 'corr22_220121aaat_data.out.corr'},
        'candname': '220121aaat',
        'output_dir': '/media/ubuntu/ssd/data/',
        'ntint': 8,
        'nfint': 1,
        'start_offset': 1716,
        'end_offset': 3252,
        'npol_out': 4}

    return params

if __name__ == '__main__':
    bin_to_uvh5()