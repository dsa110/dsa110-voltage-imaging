import numpy as np
from casacore.tables import table

DISPERSION_CONSTANT = 10/2.41

def dedisperse(vis, freq_GHz, sample_time_ms, dispersion_measure):
    """Visibilities must be (time, baseline, freq, pol)"""
    nfreqs = len(freq_GHz)
    assert vis.shape[2] == nfreqs

    dispersion_delay_ms = get_dispersion_delay_ms(freq_GHz, dispersion_measure)
    dispersion_delay_tbin = np.round(dispersion_delay_ms/sample_time_ms).astype(np.int)

    for i in range(nfreqs):
        vis[:, :, i, :] = np.roll(vis[:, :, i, :], -1*dispersion_delay_tbin[i], axis=0)

    return vis

def get_dispersion_delay_ms(freq_GHz, dispersion_measure) -> np.ndarray:
    """Calculate the dispersion delay per channel in ms."""
    reference_freq = np.max(freq_GHz)
    dispersion_delay_ms = DISPERSION_CONSTANT*dispersion_measure*(1/freq_GHz**2-1/reference_freq**2)
    return dispersion_delay_ms

def update_weights(sigma_spectrum: np.ndarray, template_filepath: str) -> None:
    """Updates a template file with weights based on pulse profile and DM.

    `sigma_spectrum` should have dimensions (Nblt, Nchans, Npols).
    """
    weight_spectrum = np.power(sigma_spectrum, -2)
    weight = weight_spectrum.sum(axis=1)
    sigma = np.power(weight, -0.5)

    with table(template_filepath) as tb:
        has_sigma_spectrum = 'SIGMA_SPECTRUM' in tb.colnames()
    if not has_sigma_spectrum:
        add_sigma_spectrum_column(template_filepath)

    with table(template_filepath, readonly=False) as tb:
        tb.putcol('SIGMA_SPECTRUM', sigma_spectrum)
        tb.putcol('WEIGHT_SPECTRUM', weight_spectrum)
        tb.putcol('WEIGHT', weight)
        tb.putcol('SIGMA', sigma)

def add_sigma_spectrum_column(template_filepath: str) -> None:
    """Add SIGMA_SPECTRUM as a column to the measurement set at `template_filepath`"""
    sigma_spectrum_desc = get_sigma_spectrum_desc()
    with table(template_filepath, readonly=False) as tb:
        tb.addcols(sigma_spectrum_desc)

def get_sigma_spectrum_desc() -> dict:
    """Generate the column description for SIGMA_SPECTRUM"""
    sigma_spectrum_desc = {'SIGMA_SPECTRUM':{
        'valueType': 'float',
        'dataManagerType': 'TiledShapeStMan',
        'dataManagerGroup': 'TiledSigmaSpectrum',
        'option': 0,
        'maxlen': 0,
        'comment': 'Estimated rms noise for each data point',
        'ndim': 2,
        '_c_order': True,
        'keywords': {}}}
    return sigma_spectrum_desc
