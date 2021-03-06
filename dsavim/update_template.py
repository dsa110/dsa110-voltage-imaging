"""Use a template to generate the measurement set from voltages."""

from typing import Tuple
from pkg_resources import resource_filename

import numpy as np
from pyuvdata import UVData
from casacore.tables import table
from astropy.time import Time
import astropy.units as u

from dsamfs.fringestopping import calc_uvw_blt
import dsautils.dsa_syslog as dsl
from dsacalib.utils import Direction

from dsavim.utils import load_params

PARAMFILE = resource_filename('dsavim', "data/voltage_corr_parameters.yaml")
PARAMS = load_params(PARAMFILE)

# Logger
LOGGER = dsl.DsaSyslogger()
LOGGER.subsystem("software")
LOGGER.app("dsacalib")

# def create_template(template_path, outms_path):
#     return
#     # We also need to include some of the things from the main table here
#     # e.g. WEIGHT_SPECTRUM, etc.
#     not_to_copy = [
#         'TIME', 'TIME_CENTROID', 'FEED TIME', 'FIELD TIME', 'OBSERVATION TIME_RANGE', 'SOURCE TIME',
#         'UVW', 'FIELD DELAY_DIR', 'FIELD PHASE_DIR', 'FIELD REFERENCE_DIR', 'SOURCE DIRECTION',
#         'DATA', 'FLAG', 'SPECTRAL_WINDOW MEAS_FREQ_REF','SPECTRAL_WINDOW CHAN_FREQ',
#         'SPECTRAL_WINDOW REF_FREQUENCY', 'SPECTRAL_WINDOW CHAN_WIDTH',
#         'SPECTRAL_WINDOW EFFECTIVE_BW', 'SPECTRAL_WINDOW RESOLUTION',
#         'SPECTRAL_WINDOW FLAG_ROW', 'SPECTRAL_WINDOW FREQ_GROUP',
#         'SPECTRAL_WINDOW FREQ_GROUP_NAME', 'SPECTRAL_WINDOW IF_CONV_CHAIN',
#         'SPECTRAL_WINDOW NAME', 'SPECTRAL_WINDOW NET_SIDEBAND',
#         'SPECTRAL_WINDOW NUM_CHAN', 'SPECTRAL_WINDOW TOTAL_BANDWIDTH', 'WEIGHT', 'SIGMA',
#         'ANTENNA1', 'ANTENNA2', asdfda]
#     with table(template) as tb:
#         default_desc = tb.getdesc()

#     items_to_change = {}
#     subtable_descs = {}
#     for key, value in default_desc['_keywords_'].items():
#         if isinstance(value, str) and value[:6] == 'Table:':
#             table_path = value.split(' ')[-1]
#             table_name = table_path.split('/')[-1]
#             with table(table_path) as tb:
#                 table_desc = tb.getdesc()

#             output_path = f'{outname}/{table_name}'
#             items_to_change[key] = f'Table: {output_path}'
#             subtable_descs[output_path] = table_desc

#     for key, value in items_to_change.items():
#         default_desc['_keywords_'][key] = value
#     del items_to_change

#     subtable_descs[outname] = default_desc
#     #with open('./default_desc.pkl', 'wb') as f:
#     #    pickle.dump(subtable_descs, f)

#     with table(outname, tabledesc=default_desc, readonly=False) as tb:
#         pass

#     for table_path, table_desc in subtable_descs.items():
#         with table(table_path, tabledesc=table_desc, readonly=False) as tb:
#             pass
#     for key in subtable_descs.keys():
#         keyname = re.findall('.ms/[_A-Z]*', key)
#         if len(keyname) == 0:
#             keyname = ''
#         else:
#             keyname = keyname[0][4:]
#         for key2 in subtable_descs[key].keys():
#             total_name = f'{keyname} {key2}'.strip()
#             if key2[0] != '_' and total_name not in not_to_copy:
#                 print(f'Update {keyname} {key2}')

def update_metadata(
        template_path: str, uvh5file: UVData, reftime_mjd: float, freq_array_Hz: np.ndarray = None,
        fringestopped: bool = False) -> None:
    """Updates a template file with real metadata.

    Questions: Should we update the first two entries of the history table?
    Currently not touching that table.
    """
    template_ms = TemplateMSMD(template_path, uvh5file.Nblts)

    if freq_array_Hz is not None:
        template_ms.update_frequency(freq_array_Hz)

    template_ms.update_obstime(convert_jd_to_mjds(uvh5file.time_array))

    ra_rad, dec_rad = get_pointing(uvh5file, reftime_mjd)
    template_ms.update_direction(ra_rad, dec_rad)

    uvw_m = calculate_uvw(uvh5file, ra_rad, dec_rad, reftime_mjd, fringestopped)
    template_ms.update_uvw(uvw_m)

def get_pointing(uvh5file: UVData, obstime_mjd: float) -> Tuple[float, float]:
    """Convert the pointing information stored in the uvh5file to J2000.

    Returns (ra, dec) in radians.
    """
    pointing = Direction(
        'HADEC', 0., uvh5file.extra_keywords['phase_center_dec'],
        obstime_mjd)
    ra_rad, dec_rad = pointing.J2000()
    return ra_rad, dec_rad

def calculate_uvw(
        uvh5file: UVData, ra_rad: float, dec_rad: float, reftime_mjd: float,
        fringestopped: bool = False) -> np.ndarray:
    """Calculate the uvw coordinates for a transit observation.

    The uvw array is calculated for the midpoint of the 2-s observation to
    decrease the time needed to create the measurement set.

    This function currently does not reproduce the uvw coordinates in the
    measurement set created through a conversion through uvfits. This is
    using the same method as used in the `uvh5_to_ms`, and may be an issue
    with precision in the uvfits file format.  Differences are on
    the order of 0.1 mm.
    """
    # TODO: Test how different these are from the values in the uvh5 files
    if fringestopped:
        print("Calculating uvws for all times")
        time_mjd = convert_jd_to_mjd(uvh5file.time_array)
        ntimes = uvh5file.Ntimes
    else:
        print(f"Calculating uvws for {reftime_mjd}")
        time_mjd = np.tile(reftime_mjd, (uvh5file.Nbls))
        ntimes = 1

    blen = calculate_blen(uvh5file)
    nblts = blen.shape[0]*ntimes

    uvw = calc_uvw_blt(
        np.tile(blen[np.newaxis, :, :], (ntimes, 1, 1)).reshape(nblts, 3),
        time_mjd,
        'J2000',
        np.tile(ra_rad*u.rad, nblts),
        np.tile(dec_rad*u.rad, nblts))

    if not fringestopped:
        uvw = np.tile(uvw[np.newaxis, :, :], (uvh5file.Ntimes, 1, 1)).reshape(uvh5file.Nblts, 3)

    # The calc_uvw_blt function returns uvw coordinates with UVData
    # sign convention.
    # We multiply this by -1 to match the CASA sign convention.
    uvw = -1.*uvw

    return uvw

def calculate_blen(uvh5file: UVData) -> np.ndarray:
    """Extract the baseline lenghts from the UVData object."""
    blen = np.zeros((uvh5file.Nbls, 3))
    for i, ant1 in enumerate(uvh5file.ant_1_array[:uvh5file.Nbls]):
        ant2 = uvh5file.ant_2_array[i]
        blen[i, ...] = uvh5file.antenna_positions[ant2, :] - uvh5file.antenna_positions[ant1, :]
    return blen

class TemplateMSMD():
    """Access and update the metatdata of the template ms."""

    def __init__(self, template_filepath: str, nrows: int):
        """Set the filepath and some other details for the template ms."""
        self.filepath = template_filepath
        self.float_type = 'float64'
        self.time_offset = -0.00011921
        self.nants = 117

        with table(self.filepath, readonly=False) as tb:
            if nrows < tb.nrows():
                rows_to_remove = tb.nrows() - nrows
                print(f"removing {rows_to_remove}")
                for i in range(rows_to_remove):
                    tb.removerows(i)
            elif nrows > tb.nrows():
                rows_to_add = nrows - tb.nrows()
                print(f"adding {rows_to_add}")
                tb.addrows(rows_to_add)

    def update_obstime(self, tobs_mjds: np.array) -> None:
        """Update the times in the template ms with the true time of observation."""
        tobs_mjds = tobs_mjds.astype(self.float_type)
        tstart = np.array(tobs_mjds[0]+self.time_offset, self.float_type)
        tstart_source = np.array(tobs_mjds[0], self.float_type)

        with table(self.filepath, readonly=False) as tb:
            tb.putcol('TIME', tobs_mjds)
            tb.putcol('TIME_CENTROID', tobs_mjds)
            tb.flush()

        with table(f"{self.filepath}/FEED", readonly=False) as tb:
            tb.putcol('TIME', np.tile(tstart, (self.nants)))
            tb.flush()

        with table(f"{self.filepath}/FIELD", readonly=False) as tb:
            tb.putcol('TIME', tstart)
            tb.flush()

        with table(f"{self.filepath}/OBSERVATION", readonly=False) as tb:
            tb.putcol('TIME_RANGE', np.tile(tstart, (1, 2)))
            tb.flush()

        with table(f"{self.filepath}/SOURCE", readonly=False) as tb:
            tb.putcol('TIME', tstart_source)
            tb.flush()

    def update_uvw(self, uvw_m: np.array) -> None:
        """Update the uvw array in the template ms with the true uvw's."""
        uvw_m = uvw_m.astype(self.float_type)

        with table(self.filepath, readonly=False) as tb:
            tb.putcol('UVW', uvw_m)
            tb.flush()

    def update_direction(self, ra_rad: float, dec_rad: float) -> None:
        """Update the directions in the template ms with the true direction."""
        direction = np.array([ra_rad, dec_rad], dtype=self.float_type)

        with table(f"{self.filepath}/FIELD", readonly=False) as tb:
            tb.putcol('DELAY_DIR', np.tile(direction, (1, 1, 1)))
            tb.putcol('PHASE_DIR', np.tile(direction, (1, 1, 1)))
            tb.putcol('REFERENCE_DIR', np.tile(direction, (1, 1, 1)))
            tb.flush()

        with table(f"{self.filepath}/SOURCE", readonly=False) as tb:
            tb.putcol('DIRECTION', np.tile(direction, (1, 1)))
            tb.flush()

    def update_frequency(self, freq_array_Hz: np.ndarray) -> None:
        """Update the frequncy information in the template ms."""
        if freq_array_Hz.ndim == 1:
            freq_array_Hz = freq_array_Hz[np.newaxis, :]
        nspw, nchan = freq_array_Hz.shape

        chan_width = np.mean(np.diff(freq_array_Hz), axis=1)

        with table(f"{self.filepath}/SPECTRAL_WINDOW", readonly=False) as tb:
            tb.putcol('MEAS_FREQ_REF', np.tile(5, nspw))
            tb.putcol('CHAN_FREQ', freq_array_Hz)
            tb.putcol('REF_FREQUENCY', freq_array_Hz[:, 0])
            tb.putcol('CHAN_WIDTH', np.tile(np.abs(chan_width), nchan))
            tb.putcol('EFFECTIVE_BW', np.tile(np.abs(chan_width), nchan))
            tb.putcol('RESOLUTION', np.tile(np.abs(chan_width), nchan))
            tb.putcol('FLAG_ROW', np.tile(False, nspw))
            tb.putcol('FREQ_GROUP', np.tile(0, nspw))
            tb.putcol('FREQ_GROUP_NAME', np.tile('none', nspw))
            tb.putcol('IF_CONV_CHAIN', np.tile(0, nspw))
            tb.putcol('NAME', np.tile('none', nspw))
            tb.putcol('NET_SIDEBAND', np.tile(1 if chan_width > 0 else 0, nspw))
            tb.putcol('NUM_CHAN', np.tile(nchan, nspw))
            tb.putcol('TOTAL_BANDWIDTH', chan_width*nchan)
            tb.flush()

class TemplateMSVis():
    """Access and update the visibilities and flags for a template ms."""

    def __init__(self, template_filepath: str, n_corr_nodes: int, vis_shape: Tuple[int, int, int]):
        """Open the template ms and instantiate the vis and flag arrays."""
        # The measurement set frequencies should be in ascending order, and there
        # should be a single spectral window.
        with table(f"{template_filepath}") as tb:
            spw = np.array(tb.DATA_DESC_ID[:])
        assert np.all(spw == spw[0])
        spw = spw[0]

        # TODO: Update to change the visibility table shape
        with table(f"{template_filepath}/SPECTRAL_WINDOW") as tb:
            freq = np.array(tb.CHAN_FREQ[:])[spw, :]
            channel_width = np.median(np.diff(freq))
            freq_ascending = channel_width > 0

        self.filepath = template_filepath
        self.vis = np.zeros(vis_shape, np.complex64)
        self.flags = np.ones(vis_shape, bool)
        self.freq = freq
        self.channel_width = channel_width
        self.freq_ascending = freq_ascending
        self.nfreq_corr = vis_shape[1]//n_corr_nodes

    def update_vis_and_flags(self, uvh5file: UVData) -> None:
        """Update the vis and flag arrays using an open uvh5 file."""
        assert uvh5file.Nfreqs == self.nfreq_corr

        uvh5_vis = uvh5file.data_array.reshape(
            uvh5file.Ntimes, uvh5file.Nbls, uvh5file.Nfreqs, uvh5file.Npols)
        uvh5_flags = uvh5file.flag_array.reshape(
            uvh5file.Ntimes, uvh5file.Nbls, uvh5file.Nfreqs, uvh5file.Npols)
        uvh5_freq = uvh5file.freq_array.squeeze(0)
        uvh5_freq_ascending = np.median(np.diff(uvh5_freq)) > 0
        print(f"uvh5_freq_ascending: {uvh5_freq_ascending}")
        print(f"template_ms.freq_ascending: {self.freq_ascending}")

        if uvh5_freq_ascending != self.freq_ascending:
            uvh5_vis = uvh5_vis[:, :, ::-1, :]
            uvh5_freq = uvh5_freq[::-1]

        # We have to use argmin here instead of searchsorted because there is a
        # 4.7e-7 Hz shift between the channels in the ms and the uvh5 file. This
        # is likely due to the numerical resolution of the fits file that is
        # currently used to convert from uvh5 to ms.
        start_chan = np.argmin(np.abs(self.freq-uvh5_freq[0]))
        end_chan = start_chan + self.nfreq_corr
        print(f"overwriting channels: {start_chan} to {end_chan}")
        print(f"uvh5 first channel: {uvh5_freq[0]}")
        print(f"template_ms corresponding channel: {self.freq[start_chan]}")
        self.vis[:, start_chan:end_chan, :] = np.conjugate(uvh5_vis.reshape(
            uvh5file.Nblts, uvh5file.Nfreqs, uvh5file.Npols))
        self.flags[:, start_chan:end_chan, :] = uvh5_flags.reshape(
            uvh5file.Nblts, uvh5file.Nfreqs, uvh5file.Npols)

    def write_vis_and_flags(self) -> None:
        """Write updated visibility and flags to the template ms."""
        with table(self.filepath, readonly=False) as tb:
            tb.putcol('DATA', self.vis)
            tb.putcol('FLAG', self.flags)
            tb.flush()

def convert_jd_to_mjds(time_jd: float) -> float:
    """Convert times between jd (Julian Date) and mjds (Modified Julian Date Seconds)."""
    time_mjd = convert_jd_to_mjd(time_jd)
    return (time_mjd*u.d).to_value(u.s)

def convert_jd_to_mjd(time_jd: float) -> float:
    """Convert times between jd (Julian Date) and mjd (Modified Julian Date)."""
    return Time(time_jd, format='jd').mjd
