"""Create uvh5 files from T3 visibilities.
"""
import os
import re
import tqdm
import h5py
import numpy as np
import astropy.units as u
import astropy.constants as c
from dsamfs.io import initialize_uvh5_file, update_uvh5_file
from dsacalib.fringestopping import calc_uvw, calc_uvw_interpolate
import dsacalib.constants as ct


class UVH5Generator:

    def __init__(
            self, candname: str, corrdir: str, declination: float, vis_params: dict, start_offset: int = None,
            end_offset: int = None):
        """
        Parameters
        ----------
        name : str
            The name of the measurement set.  If ends in 'template', won't actually read in
            any correlated data.
         pt_dec : quantity
            The pointing declination in degrees or equivalient.
        tstart : astropy.time.Time instance
            The start time of the correlated data.
        ntint : float
            The number of time bins that have been binned together (compared to the
            native correlator resolution).
        nfint : float
            The number of frequency bins to bin together before writing the ms
            (compared to the native resolution).
        filelist : dictionary
            The correlator data files for each node.
        params : dictionary
            T3 parameters.
        start_offset : ints
            The timesample to start at.  If given, end_offset must also be given.
            Defaults to transform the whole file to a ms.
        end_offset : int
        """

        self.filename = f"{corrdir}/{candname}"
        self.declination = declination*u.deg
        self.vis_params = dict(vis_params)
        self.start_offset = start_offset if start_offset is not None else 0
        self.end_offset = end_offset if end_offset is not None else len(self.vis_params['tobs'])

        # Determine which times to extract
        self.tobs = self.vis_params['tobs'][self.start_offset:self.end_offset]
        self.buvw, self.ant_bw = calculate_uvw_and_geodelay(self.vis_params, self.vis_params['tref'], self.declination)

        # Calculate parameters on the block size and frame rate for writing to the uvh5 file
        self.parse_size_parameters(self.vis_params, self.start_offset, self.end_offset)


    def process(self, corrfile: str, remove_file: bool=True) -> str:
        # Generate a uvh5 file for each corr node
        corr = re.findall(r"corr[0-9][0-9]", corrfile)[0]
        outname = f"{self.filename}_{corr}.hdf5"
        # Dont overwrite if the uvh5 file already exists
        if os.path.exists(outname):
            return outname

        fobs_corr_full = self.get_corr_frequencies(corr)
        # We only need the frequencies in the uvh5 file, so integrate
        # if we plan to integrate the data in frequency
        fobs_corr = np.median(fobs_corr_full.reshape(-1, self.vis_params['nfint']), axis=-1)

        with h5py.File(outname, 'w') as fhdf5:
            # Initialize the uvh5 file
            initialize_uvh5_file(
                fhdf5,
                self.nfreq_out,
                self.vis_params['npol'],
                self.declination.to_value(u.rad),
                self.vis_params['antenna_order'],
                self.fobs_corr,
                self.vis_params['antenna_cable_delays'])
            
            with open(corrfile, 'rb') as cfhandler:
                if self.start_offset is not None:
                    # Seek in bytes, and we have 4 bytes per item (float32)
                    cfhandler.seek(self.start_offset * 4 * self.size_params['itemspframe'])

                for i in tqdm.tqdm(range(self.size_params['nblocks'])):
                    # Define the time that we are reading in for this block
                    block_tidxs = slice(
                        i*self.size_params['framespblock'],
                        (i+1)*self.size_params['framespblock'])
                    tobs_block = self.tobs[block_tidxs]
                    buvw_block = self.buvw[block_tidxs, ...]

                    vis_chunk = self.get_visibility_chunk(cfhandler)

                    # Write the data to the uvh5 file
                    update_uvh5_file(
                        fhdf5,
                        vis_chunk.astype(np.complex64),
                        tobs_block.jd,
                        self.vis_params['tsamp'],
                        self.vis_params['bname'],
                        buvw_block,
                        np.ones(vis_chunk.shape, np.float32))

            if remove_file:
                os.remove(corrfile)

            return outname


    def calculate_uvw_and_geodelay(self, tobs: list = None, interpolate_uvws: bool = True) -> None:
        """Calculate the uvw coordinates in the correlated file.
        """
        if tobs is None:
            tobs = self.tobs

        if tobs.ndim == 0:
            interpolate_uvws = False
            ntimes = 1
        else:
            ntimes = len(tobs)

        if interpolate_uvws:
            buvw = calc_uvw_interpolate(self.vis_params['blen'], tobs, 'HADEC', 0.*u.rad, self.declination)
        else:
            bu, bv, bw = calc_uvw(
                self.vis_params['blen'], tobs.mjd, 'HADEC',
                np.tile(0.*u.rad, ntimes), np.tile(self.declination, ntimes))
            buvw = np.array([bu, bv, bw]).T

        ant_bw = buvw[:, self.vis_params['refidxs'], -1]
        buvw = np.tile(buvw, (len(tobs), 1, 1))

        return buvw, ant_bw


    def parse_size_parameters(self):
        """Define size of frames and blocks in the uvh5 file.
        """
        nchan = self.vis_params['nchan_corr']//self.vis_params['nfint']
        itemspframe = self.vis_params['nbls']*nchan*self.vis_params['npol']*2
        framespblock = 8
        if (self.end_offset - self.start_offset)%framespblock != 0:
            framespblock = self.end_offset-self.start_offset
            print(f"Changing framespblock to {framespblock}")

        itemspblock = itemspframe*framespblock
        nblocks = (self.end_offset-self.start_offset)//framespblock
        chunk_shape = (
            framespblock,
            self.vis_params['nbls'],
            nchan,
            self.vis_params['npol'])

        self.size_params = {
            'itemspframe': itemspframe,
            'framespblock': framespblock,
            'itemspblock': itemspblock,
            'nblocks': nblocks,
            'chunk_shape': chunk_shape}
        self.nfreq_out = chunk_shape[2]


    def get_corr_frequencies(self, corr: str) -> np.ndarray:
        """Determine the frequencies handled by the corr node.

        Parameters
        ----------
        corr : str
            The corr node name, e.g. 'corr00'.

        Returns
        -------
        np.ndarray
            The frequency of each channel in the visibilities.
        """
        ch0 = self.vis_params['corr_ch0'][corr]
        fobs_corr_full = self.vis_params['fobs'][ch0:(ch0+self.vis_params['nchan_corr'])]
        return fobs_corr_full


    def get_visibility_chunk(self,  cfhandler: 'FileHandler') -> np.ndarray:
        """Get a visibility chunk from the corr file open in `cfhandler`.

        Parameters
        ----------
        cfhandler : FileHandler
            The filehandler to the open correlated file.

        Returns
        -------
        np.ndarray
            The visibility chunk in complex64.
        """
        itemspblock = self.size_params['itemspblock']
        chunk_shape = self.size_params['chunk_shape']

        data = np.fromfile(
            cfhandler,
            dtype=np.float32,
            count=itemspblock
        )
        data = data.reshape(-1, 2)
        data = data[..., 0] + 1.j*data[..., 1]
        data = data.reshape(chunk_shape)

        return data


def get_total_delay(
        baseline_cable_delay: np.ndarray, ant_bw: np.ndarray, bname: list,
        antenna_order: list) -> np.ndarray:
    """Calculate total (cable plus geometric) delay for each baseline.

    Parameters
    ----------
    baseline_cable_delay : np.ndarray
        The relative cable delay on each baseline to be accounted for.
    ant_bw : np.ndarray
        The w-term for each antenna, describing the geometric delay.
    bname : list
        The name of each baseline, e.g. '24-25'.  Antenna names are indexed to 1.
    antenna_order : list
        The names of antennas, in correlator order.  Antenna names are indexed to 1.

    Returns
    -------
    np.ndarray
        The total delay for each baseline, including cable and geometric terms.
        Dimensions (time, baseline, 1, 1).
    """
    nbls = len(bname)
    ntimes = ant_bw.shape[0]
    total_delay = np.zeros((ntimes, nbls))
    for bni, bn in enumerate(bname):
        ant1, ant2 = bn.split('-')
        idx1 = antenna_order.index(int(ant1))
        idx2 = antenna_order.index(int(ant2))
        total_delay[:, bni] = baseline_cable_delay[bni] + (
            (ant_bw[:, idx1]-ant_bw[:, idx2])*u.m/c.c).to_value(u.nanosecond)
    # Reshape to match data shape
    total_delay = total_delay[:, :, np.newaxis, np.newaxis]
    return total_delay


def get_XX_YY(vis_chunk: np.ndarray) -> np.ndarray:
    """Return only the XX and YY pols of `vis_chunk`.

    Parameters
    ----------
    vis_chunk : np.ndarray
        The visbility chunk.  Last dimension is polarization, XX XY YX YY.

    Returns
    -------
    np.ndarray
        The visibility chunk with only XX YY polarizations.
    """
    return vis_chunk[..., [0, -1]]


def get_visibility_model(total_delay: np.ndarray, fobs_corr_full: np.ndarray) -> np.ndarray:
    """Construct the visibility model from the total delay.

    Parameters
    ----------
    total_delay : np.ndarray
        The total delay for each baseline, in ns.
    fobs_corr_full : np.ndarray
        The frequency of each channel in the visibilities, in GHz.

    Returns
    -------
    np.ndarray(complex64)
        The visibility model
    """
    vis_model = np.exp(2j*np.pi*fobs_corr_full[:, np.newaxis]*total_delay)
    vis_model = vis_model.astype(np.complex64)
    return vis_model


def get_mjd(armed_mjd: float, utc_start: int, specnum: int) -> float:
    """Get the start mjd of a voltage dump.

    Parameters
    ----------
    armed_mjd : float
        The time at which the snaps were armed, in mjd.
    utc_start : int
        The spectrum number at which the correlator was started.
    specnum : int
        The spectrum number of the first spectrum in the voltage dump,
        referenced to when the correlator was started.

    Returns
    -------
    tstart : float
        The start time of the voltage dump in mjd.
    """
    tstart = (armed_mjd+utc_start*4*8.192e-6/86400+
              (1/(250e6/8192/2)*specnum/ct.SECONDS_PER_DAY))
    return tstart


def calculate_uvw_and_geodelay(visparams, tobs, declination, interpolate_uvws: bool=True):

