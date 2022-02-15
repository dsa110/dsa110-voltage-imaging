import os
import shutil
import glob
import re
import numpy as np
from dsacalib.uvh5_to_ms import load_uvh5_file, set_antenna_positions, phase_visibilities
from dsacalib.uvh5_to_ms import fix_descending_missing_freqs, write_UV_to_ms
from dsaT3.dedisperse import dedisperse_and_select_times
from dsaT3.update_template import update_metadata, TemplateMSVis
from casacore.tables import tablecopy
from casatasks import virtualconcat
import astropy.units as u
from astropy.time import Time


UVH5DIR = '/media/ubuntu/ssd/data/'
CORRNAME_PATTERN = re.compile('corr[0-9][0-9]')
MSDIR = '/media/ubuntu/data/dsa110/imaging/'
TEMPLATE = f'{MSDIR}/template.ms'
VOLTAGEDIR = '/media/ubuntu/data/dsa110/T3/'
CENTRE_TIME_S = 1907*262.144e-6
DAY_TO_MS = 86400000.0
DAY_TO_S = DAY_TO_MS/1e3
REF_FREQ_GHZ = 1.530

def uvh5_to_ms(candname, candtime, dispersion_measure=None, uvh5files=None, msname=None,
               ntbins=8, centre_time_s=CENTRE_TIME_S, template_path=TEMPLATE, singlems=False):
    """Convert uvh5 to ms.

    This mostly follows the method used in the real-time system with some differences:
    * We need to account for removing geometric delays to the start of the observation.
    * We have the option to dedisperse.
    """
    if dispersion_measure is not None and template_path is not None:
        print('Currently using a template for dedispersed data is not working reliable.')
        print(f'Will create new multi-ms for {candname}.')
        template_path = None
        singlems = False

    use_template = template_path is not None

    if uvh5files is None:
        uvh5files = sorted(glob.glob(f'{UVH5DIR}{candname}_corr??.hdf5'))

    if msname is None:
        msname = f'{MSDIR}/{candname}'
    
    if os.path.exists(f'{msname}.ms'):
        shutil.rmtree(f'{msname}.ms')
    
    if (not singlems) or use_template: 
        
        if use_template:
            tablecopy(template_path, f'{msname}.ms', deep=True)
            template_ms = None
        
        for uvh5file in uvh5files:
            corr = re.findall('corr[0-9][0-9]', uvh5file)[0]
            UV, _pt_dec, ra, dec = load_uvh5_file(uvh5file)
            antenna_positions = set_antenna_positions(UV)
            process_UV(UV, dispersion_measure, ra, dec, candtime, centre_time_s, ntbins)

            if use_template:
                
                if template_ms is None:
                    update_metadata(f'{msname}.ms', UV)
                    template_ms = TemplateMSVis(
                        f'{msname}.ms', 16,
                        (UV.Nblts*UV.Nspws, UV.Nfreqs*len(uvh5files), UV.Npols))

                template_ms.update_vis_and_flags(UV)
                
            else:

                write_UV_to_ms(UV, f'{msname}_{corr}', antenna_positions)

        if use_template:

            template_ms.write_vis_and_flags()

        else:

            msnames = sorted(glob.glob(f'{msname}_corr??.ms'), reverse=True)
            virtualconcat(msnames, f'{msname}.ms')

    else:

        UV, _pt_dec, ra, dec = load_uvh5_file(uvh5files)
        antenna_positions = set_antenna_positions(UV)
        process_UV(UV, dispersion_measure, ra, dec, candtime, centre_time_s, ntbins)
        write_UV_to_ms(UV, msname, antenna_positions)

def process_UV(UV, dispersion_measure, ra, dec, candtime, centre_time_s, ntbins):
    """Phase, dedisperse and select times from UV file"""

    # TODO: account for the fact that the bws were removed for the start time
    # TODO: change ra and dec to be at the start of the burst
    # TODO: reflect that the data are actually phased in the uvh5 files

    if dispersion_measure is None:
        phase_visibilities(UV, fringestop=True, phase_ra=ra, phase_dec=dec, interpolate_uvws=True)
    else:
        phase_visibilities(UV, fringestop=False, interpolate_uvws=True)

    if dispersion_measure is not None:
        dedisperse_and_select_times_UV(UV, dispersion_measure, candtime+centre_time_s*u.s, ntbins)

    fix_descending_missing_freqs(UV)

def dedisperse_and_select_times_UV(UV, dispersion_measure, centre_time, ntbins):
    """Dedisperse the visibilities in a UVData instance.

    Dedisperse to the top channel.
    """
    freq_GHz = UV.freq_array.reshape(-1)/1e9
    nfreqs = freq_GHz.size

    tobs_jd = UV.time_array.reshape(UV.Ntimes, UV.Nbls)[:, 0]
    tobs = Time(tobs_jd, format='jd')
    centre_tbin = np.searchsorted(tobs.mjd, centre_time.mjd)
    sample_time_ms = np.mean(np.diff(tobs.mjd))*DAY_TO_MS

    vis = UV.data_array.reshape(UV.Ntimes, UV.Nbls, nfreqs, UV.Npols)
    
    vis = dedisperse_and_select_times(vis, freq_GHz, sample_time_ms, dispersion_measure, centre_tbin, ntbins, REF_FREQ_GHZ)

    UV.select(times=tobs_jd[centre_tbin-ntbins//2:centre_tbin+ntbins//2])

    UV.data_array = vis.reshape(UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols)    

def dedisperse_UV(UV, dispersion_measure):
    """Dedisperse the visibilities in a UVData instance.

    Dedisperse to the top channel.
    """
    freq_GHz = UV.freq_array.reshape(-1)/1e9
    nfreqs = freq_GHz.size

    time = UV.time_array.reshape(UV.Ntimes, UV.Nbls)[:, 0]
    sample_time_ms = np.mean(np.diff(time))*DAY_TO_MS

    vis = UV.data_array.reshape(UV.Ntimes, UV.Nbls, nfreqs, UV.Npols)
    vis = dedisperse(vis, freq_GHz, sample_time_ms, dispersion_measure, REF_FREQ_GHZ)
    UV.data_array = vis.reshape(UV.Nblts, UV.Nspws, UV.Nfreqs, UV.Npols)

def select_times_UV(UV, centre_time, ntbins):
    """Select only some times around `centre_time_s` from a UVfile.
    """
    tobs_jd = UV.time_array.reshape(UV.Ntimes, UV.Nbls)[:, 0]
    tobs = Time(tobs_jd, format='jd')
    centre_bin = np.searchsorted(tobs.mjd, centre_time.mjd)

    UV.select(times=tobs_jd[centre_bin-ntbins//2:centre_bin+ntbins//2])
