"""Use a template to generate the measurement set from voltages."""

import numpy as np
from pkg_resources import resource_filename
from pyuvdata import UVData
from casacore.tables import table
from astropy.time import Time
import astropy.units as u
import dsautils.dsa_syslog as dsl
from dsaT3.utils import load_params
from dsacalib.utils import direction as pointing_direction

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
PARAMS = load_params(PARAMFILE)

# Logger
LOGGER = dsl.DsaSyslogger()
LOGGER.subsystem("software")
LOGGER.app("dsacalib")

N_CORR_NODES = 16

def update_template(template_filepath: str, uvh5_filepaths: list):
    """Updates a template file with real data and metadata.
    
    Eventually, we will want to do this directly from the correlated data, but for
    now we use the uvh5 file as an intermediate step.  This is because it already
    properly handles outrigger delays and frequency integration, as well as reading
    data by block.
    """
    update_visibilities(template_filepath, uvh5_filepaths)
    update_metadata(template_filepath, uvh5_filepaths[0])

def update_visibilities(template_filepath: str, uvh5_filepaths: list, n_corr_nodes: int=N_CORR_NODES):
    """Updates a template file with real data.

    There are a few ways that we can do this:
    The most error-proof way is to use casacore, but it may not be able to
    handle large data sizes.

    We can stream bytes directly from the input file to the ms, as long as
    we get the offsets correct.

    A faster way may be to write the file directly to the visibility table,
    although this depends on the speed of writing the correlator output to
    the spinning disk (where the measurement sets are stored) and the ssd
    (where the voltages are stored and correlated).

    For now we will use casacore.
    """
    template_ms = TemplateMSVis(template_filepath, n_corr_nodes)

    for uvh5_filepath in uvh5_filepaths:
        UV = UVData()
        UV.read(uvh5_filepath, file_type='uvh5')
        template_ms.update_vis_and_flags(UV)

    template_ms.write_vis_and_flags()

def update_metadata(template_path: str, uvh5_filepath: str) -> None:
    """Updates a template file with real metadata.

    Questions: Should we update the first two entries of the history table?
    Currently not touching that table.
    """
    template_ms = TemplateMSMD(template_path)

    UV = UVData()
    UV.read(uvh5_filepath, file_type='uvh5')

    template_ms.update_obstime(convert_jd_to_mjds(UV.time_array))

    template_ms.update_uvw(UV.uvw_array)

    ra_rad, dec_rad = get_pointing(UV)
    template_ms.update_direction(ra_rad, dec_rad)

def get_pointing(uvh5file: UVData) -> tuple:
    """Convert the pointing information stored in the uvh5file to J2000."""
    obstime_mjd = Time(np.mean(uvh5file.time_array), format='jd').mjd
    pointing = pointing_direction(
        'HADEC', 0., uvh5file.extra_keywords['phase_center_dec'],
        obstime_mjd)
    ra_rad, dec_rad = pointing.J2000()
    return ra_rad, dec_rad
    
def update_weights(simulated_dynspec: np.ndarray, template_filepath: str) -> None:
    """Updates a template file with weights based on pulse profile and DM."""
    pass

class TemplateMSMD():
    """Access and update the metatdata of the template ms."""

    def __init__(self, template_filepath: str):
        """Set the filepath and some other details for the template ms."""
        self.filepath = template_filepath
        self.float_type = 'float64'
        self.time_offset = -0.00011921
        self.nants = 117

    def update_obstime(self, tobs_mjds: np.array):
        """Update the times in the template ms with the true time of observation."""
        tobs_mjds = tobs_mjds.astype(self.float_type)
        tstart = np.array(tobs_mjds[0]+self.time_offset, self.float_type)
        tstart_source = np.array(tobs_mjds[0], self.float_type)

        with table(self.filepath, readonly=False) as tb:
            tb.putcol('TIME', tobs_mjds)
            tb.putcol('TIME_CENTROID', tobs_mjds)

        with table(f'{self.filepath}/FEED', readonly=False) as tb:
            tb.putcol('TIME', np.tile(tstart, (self.nants)))

        with table(f'{self.filepath}/FIELD', readonly=False) as tb:
            tb.putcol('TIME', tstart)

        with table(f'{self.filepath}/OBSERVATION', readonly=False) as tb:
            tb.putcol('TIME_RANGE', np.tile(tstart, (1, 2)))

        with table(f'{self.filepath}/SOURCE', readonly=False) as tb:
            tb.putcol('TIME', tstart_source)

    def update_uvw(self, uvw_m: np.array):
        """Update the uvw array in the template ms with the true uvw's."""
        uvw_m = uvw_m.astype(self.float_type)

        with table(self.filepath, readonly=False) as tb:
            tb.putcol('UVW', -1*uvw_m)

    def update_direction(self, ra_rad: float, dec_rad: float):
        """Update the directions int he template ms with the true direction."""
        direction = np.array([ra_rad, dec_rad], dtype=self.float_type)

        with table(f'{self.filepath}/FIELD', readonly=False) as tb:
            tb.putcol('DELAY_DIR', np.tile(direction, (1, 1, 1)))
            tb.putcol('PHASE_DIR', np.tile(direction, (1, 1, 1)))
            tb.putcol('REFERENCE_DIR', np.tile(direction, (1, 1, 1)))

        with table(f'{self.filepath}/SOURCE', readonly=False) as tb:
            tb.putcol('DIRECTION', np.tile(direction, (1, 1)))

class TemplateMSVis():
    """Access and update the visibilities and flags for a template ms."""

    def __init__(self, template_filepath: str, n_corr_nodes: int=N_CORR_NODES):
        """Open the template ms and instantiate the vis and flag arrays."""
        # The measurement set frequencies should be in ascending order, and there
        # should be a single spectral window.
        with table(f'{template_filepath}/SPECTRAL_WINDOW') as tb:
            freq = np.array(tb.CHAN_FREQ[:]).squeeze(0)
            assert np.median(np.diff(freq)) > 0

        with table(template_filepath) as tb:
            vis = np.array(tb.DATA[:])
            vis_shape = vis.shape
            vis_dtype = vis.dtype

        self.filepath = template_filepath
        self.vis = np.zeros(vis_shape, vis_dtype)
        self.flags = np.ones(vis.shape, bool)
        self.freq = freq
        self.nfreq_corr = vis.shape[1]//n_corr_nodes

    def update_vis_and_flags(self, uvh5file: UVData):
        """Update the vis and flag arrays using an open uvh5 file."""
        assert uvh5file.Nfreqs == self.nfreq_corr

        uvh5_vis = uvh5file.data_array.reshape(uvh5file.Ntimes, uvh5file.Nbls, uvh5file.Nfreqs, uvh5file.Npols)
        uvh5_flags = uvh5file.flag_array.reshape(uvh5file.Ntimes, uvh5file.Nbls, uvh5file.Nfreqs, uvh5file.Npols)
        uvh5_freq = uvh5file.freq_array.squeeze(0)

        if np.median(np.diff(uvh5_freq)) < 0:
            uvh5_vis = uvh5_vis[:, :, ::-1, :]
            uvh5_freq = uvh5_freq[::-1]

        # We have to use argmin here instead of searchsorted because there is a
        # 4.7e-7 Hz shift between the channels in the ms and the uvh5 file. This
        # is likely due to the numerical resolution of the fits file that is
        # currently used to convert from uvh5 to ms.
        start_chan = np.argmin(np.abs(self.freq-uvh5_freq[0]))
        end_chan = start_chan + self.nfreq_corr
        self.vis[:, start_chan:end_chan, :] = np.conjugate(uvh5_vis.reshape(uvh5file.Nblts, uvh5file.Nfreqs, uvh5file.Npols))
        self.flags[:, start_chan:end_chan, :] = uvh5_flags.reshape(uvh5file.Nblts, uvh5file.Nfreqs, uvh5file.Npols)

    def write_vis_and_flags(self):
        """Write updated visibility and flags to the template ms."""
        with table(self.filepath, readonly=False) as tb:
            tb.putcol('DATA', self.vis)
            tb.putcol('FLAG', self.flags)

def convert_jd_to_mjds(time_jd):
    """Convert times between jd (Julian Date) and mjds (Modified Julian Date Seconds)."""
    return (Time(time_jd, format='jd').mjd*u.d).to_value(u.s)
