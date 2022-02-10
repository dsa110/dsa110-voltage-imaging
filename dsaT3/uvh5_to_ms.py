import glob
import re
import os
import shutil
import numpy as np
from dsacalib.uvh5_to_ms import load_uvh5_file, set_antenna_positions, phase_visibilities
from dsaT3.dedisperse import add_sigma_spectrum_column
from dsaT3.imaging import calibrate_T3ms_percorrnode
from dsaT3.utils import find_beamformer_weights
from casacore.tables import table, addImagingColumns
from casatasks import virtualconcat

UVH5DIR = '/media/ubuntu/ssd/data/'
CORRNAME_PATTERN = re.compile('corr[0-9][0-9]')
MSDIR = '/media/ubuntu/data/dsa110/imaging/'
VOLTAGEDIR = '/media/ubuntu/data/dsa110/T3/'

def uvh5_to_ms(candname, candtime, uvh5files=None):
    if uvh5files is None:
        uvh5files = sorted(glob.glob(f'{UVH5DIR}{candname}_corr??.hdf5'))

    for file in uvh5files:
        print(f'processing {file}')
        filename = file.split('/')[-1]
        corrname = CORRNAME_PATTERN.findall(filename)[0]
        mspath = f'{MSDIR}{candname}_{corrname}.ms'
        if os.path.exists(mspath):
            shutil.rmtree(mspath)

        UV, _pt_dec, ra, dec = load_uvh5_file(file)
        antenna_positions = set_antenna_positions(UV)
        phase_visibilities(UV, fringestop=True, phase_ra=ra, phase_dec=dec, interpolate_uvws=True)

        print(f'writing {file}')
        # Swap uvw coordinates and visibilities
        # TODO: Investigate why pyuvdata doesn't do this automatically
        UV.uvw_array = -1*UV.uvw_array
        UV.data_array = np.conjugate(UV.data_array)

        UV.write_ms(str(mspath))

        with table(f'{mspath}/ANTENNA', readonly=False) as tb:
            tb.putcol('POSITION', antenna_positions)

        addImagingColumns(mspath)

        with table(mspath) as tb:
            has_sigma_spectrum = 'SIGMA_SPECTRUM' in tb.colnames()
        if not has_sigma_spectrum:
            add_sigma_spectrum_column(mspath)

    msnames = {CORRNAME_PATTERN.findall(msname)[0]: str(msname)[:-3]
               for msname in glob.glob(f'{MSDIR}{candname}_corr??.ms')}

    print('calibrating')
    bfweights = find_beamformer_weights(candtime)
    calibrate_T3ms_percorrnode(msnames, bfweights)

    print('concatenating')
    virtualconcat(sorted([f'{msname}.ms' for msname in msnames.values()]), concatvis=f'{MSDIR}{candname}.ms')
