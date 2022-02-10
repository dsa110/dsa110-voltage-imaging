"""Write to a template ms from T3 visibilities.
"""
import os
import re
import tqdm
import numpy as np
from pkg_resources import resource_filename
import astropy.units as u
import astropy.constants as c
from antpos.utils import get_itrf
from dsamfs.io import initialize_uvh5_file, update_uvh5_file
from dsacalib.fringestopping import calc_uvw
import dsacalib.constants as ct
from dsaT3.utils import load_params
from dsaT3.generate_uvh5 import parse_visibility_parameters, parse_size_parameters, calculate_uvw_and_geodelay
from dsaT3.generate_uvh5 import get_total_delay, get_visibility_chunk, get_XX_YY

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
T3PARAMS = load_params(PARAMFILE)
NPOL_OUT = 2

def corr_to_ms(corrfile, measurement_set, pt_dec, tstart, ntint, nfint, params=T3PARAMS,
               start_offset=0, end_offset=-1, npolout=NPOL_OUT):
    """Converts a binary correlated file to a measurement set.
    
    The measurement set should already be created.
    """
    if end_offset < 0:
        end_offset = params['nsubint']//ntint

    # Parse further parameters
    vis_params = parse_visibility_parameters(params, tstart, ntint)

    # Determine which times to extract
    tobs = vis_params['tobs'][start_offset:end_offset]

    # Calculate parameters on the block size and frame rate for writing to the uvh5 file
    size_params = parse_size_parameters(vis_params, start_offset, end_offset, nfint, npol_out)
    nfreq_out = size_params['output_chunk_shape'][2]

    # Calculate the frequencies
    corr = re.findall('corr[0-9][0-9]', corrfile)[0]
    fobs_corr_full = get_corr_frequencies(vis_params, corr)
    fobs_corr = np.median(fobs_corr_full.reshape(-1, nfint), axis=-1)
    # lamb = None
    
    
    # Get the meridian uvw coordinates, which are what we will write out
    # blen, bname = get_blen()
    # uvw_m = None
    
    # The measurement set should already be created
    # Initialize the chunk size, times, frequencies and polarizations
    with open(corrfile, 'rb') as cfhandler:
        # Seek in the corr file
        cfhandler.seek(start_offset*4*size_params['itemspframe'])
        for i in tqdm.tqdm(size_params['nblocks']):
            tobs_block = tobs[
                i*size_params['framespblock']:(i+1)*size_params['framespblock']]

            buvw, ant_bw = calculate_uvw_and_geodelay(
                vis_params,
                size_params,
                tobs_block, pt_dec,
                interpolate_uvws=True)
            total_delay = get_total_delay(
                vis_params['baseline_cable_delays'],
                ant_bw,
                vis_params['bname'],
                vis_params['antenna_order'])

            # Read a block of data from the correlated file.
            vis_chunk = get_visibility_chunk(
                cfhandler,
                size_params['itemspblock'],
                size_params['input_chunk_shape'])
            if vis_params['npol'] == 4 and npol_out == 2:
                vis_chunk = get_XX_YY(vis_chunk)

            # Apply outrigger and geometric delays
            # We won't perform any phasing here - we just write the data
            # directly to the uvh5 file.
            vis_model = get_visibility_model(total_delay, fobs_corr_full)
            vis_chunk /= vis_model

            # Squeeze the data in frequency - this has to be done *after* applying
            # the delays
            if nfint > 1:
                vis_chunk = vis_chunk.reshape(
                    size_params['output_chunk_shape'][0],
                    size_params['output_chunk_shape'][1],
                    size_params['output_chunk_shape'][2],
                    nfint,
                    size_params['output_chunk_shape'][3]
                ).mean(axis=3)

            # Phase to the meridian
            # phase_model = generate_phase_model_antbased(
            #     uvw_m, buvws, vis_params['nbls'], size_params['framespblock'], lamb, UV.ant_1_array[:UV.Nbls],
            #     UV.ant_2_array[:UV.Nbls])
            # vis_chunk /= phase_model

            # Update the appropriate spot of the measurement set
            # TODO: chunk approach to measurement set updating.
        