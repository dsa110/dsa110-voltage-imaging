import glob
import re
import numpy as np
from dsacalib.uvh5_to_ms import load_uvh5_file, set_antenna_positions, phase_visibilities
from dsacalib.uvh5_to_ms import fix_descending_missing_freqs, write_UV_to_ms
from dsaT3.dedisperse import dedisperse

UVH5DIR = '/media/ubuntu/ssd/data/'
CORRNAME_PATTERN = re.compile('corr[0-9][0-9]')
MSDIR = '/media/ubuntu/data/dsa110/imaging/'
VOLTAGEDIR = '/media/ubuntu/data/dsa110/T3/'
CENTRE_TIME_S = 1907*262.144e-6
DAY_TO_MS = 86400000.0
DAY_TO_S = DAY_TO_MS/1e3

def uvh5_to_ms(candname, dispersion_measure=None, uvh5files=None, msname=None,
               ntbins=128, centre_time_s=CENTRE_TIME_S):
    """Convert uvh5 to ms.

    This mostly follows the method used in the real-time system with some differences:
    * We need to account for removing geometric delays to the start of the observation.
    * We have the option to dedisperse.
    """
    if uvh5files is None:
        uvh5files = sorted(glob.glob(f'{UVH5DIR}{candname}_corr??.hdf5'))

    if msname is None:
        msname = f'{MSDIR}/{candname}.ms'

    UV, _pt_dec, ra, dec = load_uvh5_file(uvh5files)
    antenna_positions = set_antenna_positions(UV)

    # TODO: account for the fact that the bws were removed for the start time
    phase_visibilities(UV, fringestop=True, phase_ra=ra, phase_dec=dec, interpolate_uvws=True)

    if dispersion_measure is not None:
        dedisperse_UV(UV, dispersion_measure)
        select_times_UV(UV, centre_time_s, ntbins)

    fix_descending_missing_freqs(UV)
    write_UV_to_ms(UV, msname, antenna_positions)

def dedisperse_UV(UV, dispersion_measure):
    """Dedisperse the visibilities in a UVData instance.

    Dedisperse to the top channel.
    """
    freq_GHz = UV.freq_array.reshape(-1)/1e9
    nfreqs = freq_GHz.size

    time = UV.time_array.reshape(UV.Ntimes, UV.Nbls)[:, 0]
    sample_time_ms = np.mean(np.diff(time))*DAY_TO_MS

    vis = UV.data_array.reshape(UV.Ntimes, UV.Nbls, nfreqs, UV.Npols)
    vis = dedisperse(vis, freq_GHz, sample_time_ms, dispersion_measure)
    UV.data_array = vis.reshape(UV.Ntimes, UV.Nbls, UV.Nspws, UV.Nfreqs, UV.Npols)

def select_times_UV(UV, centre_time_s, ntbins):
    """Select only some times around `centre_time_s` from a UVfile.
    """
    tobs = UV.time_array.reshape(UV.Ntimes, UV.Nbls)[:, 0]
    tobs_s = (tobs - tobs[0])*DAY_TO_S
    centre_bin = np.searchsorted(tobs_s, centre_time_s)

    UV.select(times=tobs[centre_bin-ntbins//2:centre_bin+ntbins//2])

# def uvh5_to_ms(candname, candtime, uvh5files=None):
#     if uvh5files is None:
#         uvh5files = sorted(glob.glob(f'{UVH5DIR}{candname}_corr??.hdf5'))

#     for file in uvh5files:
#         print(f'processing {file}')
#         filename = file.split('/')[-1]
#         corrname = CORRNAME_PATTERN.findall(filename)[0]
#         mspath = f'{MSDIR}{candname}_{corrname}.ms'
#         if os.path.exists(mspath):
#             shutil.rmtree(mspath)

#         UV, _pt_dec, ra, dec = load_uvh5_file(file)
#         antenna_positions = set_antenna_positions(UV)
#         phase_visibilities(UV, fringestop=True, phase_ra=ra, phase_dec=dec, interpolate_uvws=True)

#         print(f'writing {file}')
#         # Swap uvw coordinates and visibilities
#         # TODO: Investigate why pyuvdata doesn't do this automatically
#         UV.uvw_array = -1*UV.uvw_array
#         UV.data_array = np.conjugate(UV.data_array)

#         UV.write_ms(str(mspath))

#         with table(f'{mspath}/ANTENNA', readonly=False) as tb:
#             tb.putcol('POSITION', antenna_positions)

#         addImagingColumns(mspath)

#         with table(mspath) as tb:
#             has_sigma_spectrum = 'SIGMA_SPECTRUM' in tb.colnames()
#         if not has_sigma_spectrum:
#             add_sigma_spectrum_column(mspath)

#     msnames = {CORRNAME_PATTERN.findall(msname)[0]: str(msname)[:-3]
#                for msname in glob.glob(f'{MSDIR}{candname}_corr??.ms')}

#     print('calibrating')
#     bfweights = find_beamformer_weights(candtime)
#     calibrate_T3ms_percorrnode(msnames, bfweights)

#     print('concatenating')
#     virtualconcat(
#         sorted([f'{msname}.ms' for msname in msnames.values()]),
#         concatvis=f'{MSDIR}{candname}.ms')
