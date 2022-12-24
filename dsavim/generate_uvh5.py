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
from antpos import utils

def generate_uvh5(
        name: str, pt_dec: 'astropy.Quantity', corrfile: str, vis_params: dict,
        start_offset: int = None, end_offset: int = None) -> str:
    """Generates a measurement set from the T3 correlations.

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
        The timesample to end at.

    Returns
    -------
    str
        The name of the measurement set created.
    """

    if start_offset is None:
        start_offset = 0
    if end_offset is None:
        end_offset = len(vis_params['tobs'])

    # Determine which times to extract
    tobs = vis_params['tobs'][start_offset:end_offset]
    buvw, _ant_bw = calculate_uvw_and_geodelay(vis_params, vis_params['tref'], pt_dec)
    buvw = np.tile(buvw, (len(tobs), 1, 1))

    # Calculate parameters on the block size and frame rate for writing to the uvh5 file
    size_params = parse_size_parameters(vis_params, start_offset, end_offset)
    nfreq_out = size_params['chunk_shape'][2]

    # Generate a uvh5 file for each corr node
    sb = re.findall(r"sb[0-9][0-9]", corrfile)[0]
    _corrlist = {'sb00':'corr03','sb01':'corr04','sb02':'corr05','sb03':'corr06','sb04':'corr07',
                 'sb05':'corr08','sb06':'corr10','sb07':'corr11','sb08':'corr12','sb09':'corr14',
                 'sb10':'corr15','sb11':'corr16','sb12':'corr18','sb13':'corr19','sb14':'corr21',
                 'sb15':'corr22'}
    corr = _corrlist[sb]
    
    outname = f"{name}_{sb}.hdf5"
    # Dont overwrite if the uvh5 file already exists
    if os.path.exists(outname):
        return outname

    fobs_corr_full = get_corr_frequencies(vis_params, corr)
    # We only need the frequencies in the uvh5 file, so integrate
    # if we plan to integrate the data in frequency
    fobs_corr = np.median(fobs_corr_full.reshape(-1, vis_params['nfint']), axis=-1)

    with h5py.File(outname, 'w') as fhdf5:
        # Initialize the uvh5 file
        initialize_uvh5_file(
            fhdf5,
            nfreq_out,
            vis_params['npol'],
            pt_dec.to_value(u.rad),
            vis_params['antenna_order'],
            fobs_corr,
            vis_params['snapdelays'],
            vis_params['ant_itrf'],
            vis_params['nants_telescope'])
        with open(corrfile, 'rb') as cfhandler:
            if start_offset is not None:
                # Seek in bytes, and we have 4 bytes per item (float32)
                cfhandler.seek(start_offset*4*size_params['itemspframe'])

            for i in tqdm.tqdm(range(size_params['nblocks'])):
                # Define the time that we are reading in for this block
                block_tidxs = slice(
                    i*size_params['framespblock'],
                    (i+1)*size_params['framespblock'])
                tobs_block = tobs[block_tidxs]
                buvw_block = buvw[block_tidxs, ...]

                vis_chunk = get_visibility_chunk(
                    cfhandler,
                    size_params['itemspblock'],
                    size_params['chunk_shape'])

                # Write the data to the uvh5 file
                update_uvh5_file(
                    fhdf5,
                    vis_chunk.astype(np.complex64),
                    tobs_block.jd,
                    vis_params['tsamp'],
                    vis_params['bname'],
                    buvw_block,
                    np.ones(vis_chunk.shape, np.float32))

    return outname

def parse_size_parameters(vis_params: dict, start_offset: int, end_offset: int) -> dict:
    """Define size of frames and blocks in the uvh5 file.

    Parameters
    ----------
    vis_params : dict
        Parameters that describe the visibilities to be read in.
    start_offset, end_offset : int
        Time bins start_offset:end_offset will be read from the visibilities.
    nfint : int
        The number of frequency bins to integrate by before writing the visibilities.

    Returns
    -------
    size_params : dict
        Parameters that describe the size of data chunks to be read in and written out.
    """
    nchan = vis_params['nchan_corr']//vis_params['nfint']
    itemspframe = vis_params['nbls']*nchan*vis_params['npol']*2
    framespblock = 8
    if (end_offset - start_offset)%framespblock != 0:
        framespblock = end_offset-start_offset
        print(f"Changing framespblock to {framespblock}")

    itemspblock = itemspframe*framespblock
    nblocks = (end_offset-start_offset)//framespblock
    chunk_shape = (
        framespblock,
        vis_params['nbls'],
        nchan,
        vis_params['npol'])
    size_params = {
        'itemspframe': itemspframe,
        'framespblock': framespblock,
        'itemspblock': itemspblock,
        'nblocks': nblocks,
        'chunk_shape': chunk_shape}
    return size_params

def calculate_uvw_and_geodelay(
        vis_params: dict, tobs: np.ndarray, pt_dec: float, interpolate_uvws: bool = True) -> tuple:
    """Calculate the uvw coordinates in the correlated file.

    Parameters
    ----------
    vis_params : dict
        Parameters describing the correlated visibilities.
    size_params : dict
        Parameters describing the size of the data to be written to uvh5.
    tobs : np.ndarray(float)
        The time of each time bin in the visbilities in mjd.
    pt_dec : float
        The pointing declination of the observation in rad.

    Returns
    -------
    buvw : np.ndarray(float)
        The uvw coordinates of each baseline, dimensions (time, baseline, 3).
    ant_bw : np.ndarray(float)
        The geometric path length, w, to each antenna, dimensions (time, antenna).
    """
    if tobs.ndim == 0:
        interpolate_uvws = False
        ntimes = 1
    else:
        ntimes = len(tobs)

    if interpolate_uvws:
        buvw = calc_uvw_interpolate(vis_params['blen'], tobs, 'HADEC', 0.*u.rad, pt_dec)
    else:
        bu, bv, bw = calc_uvw(
            vis_params['blen'], tobs.mjd, 'HADEC',
            np.tile(0.*u.rad, ntimes), np.tile(pt_dec, ntimes))
        buvw = np.array([bu, bv, bw]).T

    ant_bw = buvw[:, vis_params['refidxs'], -1]

    return buvw, ant_bw

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
        total_delay[:, bni] = 1.0*baseline_cable_delay[bni] + (
            (ant_bw[:, idx1]-ant_bw[:, idx2])*u.m/c.c).to_value(u.nanosecond)
    # Reshape to match data shape
    total_delay = total_delay[:, :, np.newaxis, np.newaxis]
    return total_delay

def get_corr_frequencies(vis_params: dict, corr: str) -> np.ndarray:
    """Determine the frequencies handled by the corr node.

    Parameters
    ----------
    vis_params : dict
        Parameters describing the correlated visibilities to be read in.
    corr : str
        The corr node name, e.g. 'corr00'.

    Returns
    -------
    np.ndarray
        The frequency of each channel in the visibilities.
    """

    ch0 = vis_params['corr_ch0'][corr]
    fobs_corr_full = vis_params['fobs'][ch0:(ch0+vis_params['nchan_corr'])]
    return fobs_corr_full

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

def get_visibility_chunk(
        cfhandler: 'FileHandler', itemspblock: int, chunk_shape: tuple) -> np.ndarray:
    """Get a visibility chunk from the corr file open in `cfhandler`.

    Parameters
    ----------
    cfhandler : FileHandler
        The filehandler to the open correlated file.
    itemspblock : int
        The number of float32's per block in the correlated file.
    chunk_shape : tuple
        The desired shape of the chunk read in.

    Returns
    -------
    np.ndarray
        The visibility chunk in complex64.
    """
    data = np.fromfile(
        cfhandler,
        dtype=np.float32,
        count=itemspblock
    )
    data = data.reshape(-1, 2)
    data = data[..., 0] + 1.j*data[..., 1]
    data = data.reshape(chunk_shape)
    return data

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
