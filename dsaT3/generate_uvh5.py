"""Create uvh5 files from T3 visibilities.
"""
import os
import h5py
import numpy as np
from pkg_resources import resource_filename
import astropy.units as u
import astropy.constants as c
from antpos.utils import get_itrf
from dsamfs.io import initialize_uvh5_file, update_uvh5_file
from dsacalib.fringestopping import calc_uvw
import dsacalib.constants as ct
from dsaT3.utils import load_params

PARAMFILE = resource_filename('dsaT3', 'data/T3_parameters.yaml')
T3PARAMS = load_params(PARAMFILE)
NPOL_OUT = 2

def generate_uvh5(name, pt_dec, tstart, ntint, nfint, filelist, params=T3PARAMS,
                  start_offset=None, end_offset=None) -> str:
    """Generates a measurement set from the T3 correlations.

    Parameters
    ----------
    name : str
        The name of the measurement set.
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
        end_offset = params['nsubint']//ntint

    # Parse further parameters
    vis_params = parse_visibility_parameters(params, tstart, ntint)

    # Determine which times to extract
    tobs = vis_params['tobs'][start_offset:end_offset]

    # Calculate parameters on the block size and frame rate for writing to the uvh5 file
    size_params = parse_size_parameters(vis_params, start_offset, end_offset, nfint)
    nfreq_out = size_params['output_chunk_shape'][2]

    # Generate a uvh5 file for each corr node
    for corr, corrfile in filelist.items():

        outname = '{1}_{0}.hdf5'.format(corr, name)
        # Dont overwrite if the uvh5 file already exists
        if os.path.exists(outname):
            break

        fobs_corr_full = get_corr_frequencies(vis_params, corr)
        # We only need the frequencies in the uvh5 file, so integrate
        # if we plan to integrate the data in frequency
        fobs_corr = np.median(fobs_corr_full.reshape(-1, nfint), axis=-1)

        with h5py.File(outname, 'w') as fhdf5:
            # Initialize the uvh5 file
            initialize_uvh5_file(
                fhdf5,
                nfreq_out,
                NPOL_OUT,
                pt_dec.to_value(u.rad),
                vis_params['antenna_order'],
                fobs_corr,
                vis_params['antenna_cable_delays']
            )

            # Open the correlated file and seek the desired start position
            # If we're making a template we can't open this file
            # How can we open it with "with" if we want the alternative to do nothing?
            # Opening it each time seems like too much overhead
            with open(corrfile, 'rb') as cfhandler:
                if start_offset is not None:
                    cfhandler.seek(start_offset*32*size_params['itemspframe'])

                for i in range(size_params['nblocks']):
                    # Define the time that we are reading in for this block
                    tobs_block = tobs[
                        i*size_params['framespblock']:(i+1)*size_params['framespblock']]

                    buvw, ant_bw = calculate_uvw_and_geodelay(
                        vis_params,
                        size_params,
                        tobs_block, pt_dec)
                    total_delay = get_total_delay(
                        vis_params['baseline_cable_delay'],
                        ant_bw,
                        vis_params['bname'],
                        vis_params['antenna_order'])

                    # If we're making a template we can't read the data
                    # Read a block of data from the correlated file.
                    vis_chunk = get_visibility_chunk(
                        cfhandler,
                        size_params['itemspblock'],
                        size_params['input_chunk_shape'])
                    if vis_params['npols'] == 4 and NPOL_OUT == 2:
                        vis_chunk = vis_chunk.get_XX_YY(vis_chunk)

                    # Apply outrigger and geometric delays
                    # We won't perform any phasing here - we just write the data
                    # directly to the uvh5 file.
                    vis_model = get_visibility_model(total_delay, fobs_corr_full)
                    data /= vis_model
                    # Squeeze the data in frequency - this has to be done *after* applying
                    # the delays
                    if nfint > 1:
                        data = data.reshape(
                            size_params['output_chunk_shape'][0],
                            size_params['output_chunk_shape'][1],
                            size_params['output_chunk_shape'][2],
                            nfint,
                            size_params['output_chunk_shape'][3]
                        ).mean(axis=3)

                    # If we are making a template, we want this data instead:
                    # data = np.zeros((framespblock, nbls, len(fobs_corr), nfint, 2), dtype=np.complex64)

                    # Write the data to the uvh5 file
                    update_uvh5_file(
                        fhdf5,
                        data.astype(np.complex64),
                        tobs_block.jd,
                        vis_params['tsamp'],
                        vis_params['bname'],
                        buvw,
                        np.ones(data.shape, np.float32)
                    )

    return outname

def parse_size_parameters(vis_params: dict, start_offset: int, end_offset: int,
                          nfint: int) -> dict:
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
    itemspframe = vis_params['nbls']*vis_params['nchan_corr']*vis_params['npol']*2
    framespblock = 2
    itemspblock = itemspframe*framespblock
    assert (end_offset - start_offset)%framespblock == 0
    nblocks = (end_offset-start_offset)//framespblock
    input_chunk_shape = (
        framespblock,
        vis_params['nbls'],
        vis_params['nchan_corr'],
        vis_params['npol'])
    output_chunk_shape = (framespblock, vis_params['nbls'], vis_params['nfreq']//nfint, NPOL_OUT)
    size_params = {
        'framespblock': framespblock,
        'itemspblock': itemspblock,
        'nblocks': nblocks,
        'input_chunk_shape': input_chunk_shape,
        'output_chunk_shape': output_chunk_shape}
    return size_params

def parse_visibility_parameters(params: dict, tstart: 'astropy.time.Time', ntint: int) -> dict:
    """Parse details about the data to be extracted.

    Parameters
    ----------
    params : dict
        Parameters passed by the user definining the correlator and T3 system.
    tstart : astropy.time.Time
        Start time of the correlated visibility file.
    ntint : int
        The number of time bins integrated together in the correlated data
        compared to the native resolution.

    Returns
    -------
    dict
        Parameters that describe the visibilities to be read in.
    """
    antenna_order = params['antennas']
    fobs = params['f0_GHz']+params['deltaf_MHz']*1e-3*np.arange(params['nchan'])
    nant = len(antenna_order)
    nbls = (nant*(nant+1))//2

    # Visibilities have already been integrated in time, so account for that when
    # determining the observation times.
    tsamp = params['deltat_s']*ntint*u.s
    tobs = tstart + (np.arange(params['nsubint']//ntint)+0.5)*tsamp

    # Get baselines
    blen, bname = get_blen(antenna_order)
    # Get indices for baselines to the reference antenna
    refidxs = []
    refant = str(antenna_order[0])
    for i, bn in enumerate(bname):
        if refant in bn.split('-'):
            refidxs += [i]

    cable_delays = get_cable_delays(params['outrigger_delays'], bname)

    vis_params = {
        # Baselines
        'antenna_order': antenna_order,
        'blen': blen,
        'bname': bname,
        'nbls': nbls,
        'refidxs': refidxs,
        'antenna_cable_delays': params['outrigger_delays'],
        'baseline_cable_delays': cable_delays,
        # Time
        'tsamp': tsamp,
        'tobs': tobs,
        # Frequency
        'fobs': fobs,
        'corr_ch0': params['ch0'],
        'nchan_corr': params['nchan_corr'],
        # Polarization
        'npol': params['npol']}

    return vis_params

def get_cable_delays(outrigger_delays: dict, bname: list) -> np.ndarray:
    """Calculate cable delays from the measured outrigger cable delays.

    Antenna names in both input parameters are indexed at 1.

    Parameters
    ----------
    outrigger_delays : dict
        The delay for each antenna.  Missing keys have 0 delay.
    bname : list
        The name of each baseline.

    Returns
    -------
    np.ndarray
        The cable delay for each baseline in bname.
    """
    delays = np.zeros(len(bname), dtype=np.int)
    for i, bn in enumerate(bname):
        ant1, ant2 = bn.split('-')
        delays[i] = outrigger_delays.get(int(ant1), 0)-\
                    outrigger_delays.get(int(ant2), 0)
    return delays

def calculate_uvw_and_geodelay(vis_params: dict, size_params: dict,
                               tobs: np.ndarray, pt_dec: float) -> tuple:
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
    bu, bv, bw = calc_uvw(
        vis_params['blen'],
        tobs.mjd,
        'HADEC',
        np.zeros(size_params['framespblock'])*u.rad,
        np.ones(size_params['framespblock'])*pt_dec
    )
    buvw = np.array([bu, bv, bw]).T
    ant_bw = bw[vis_params['refidxs']].T
    return buvw, ant_bw

def get_total_delay(baseline_cable_delay: np.ndarray, ant_bw: np.ndarray, bname: list,
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

def get_visibility_chunk(cfhandler: "FileHandler", itemspblock: int,
                         chunk_shape: tuple) -> np.ndarray:
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

def get_blen(antennas: list) -> tuple:
    """Gets the baseline lengths for a subset of antennas.

    Parameters
    ----------
    antennas : list
        The antennas used in the array.

    Returns
    -------
    blen : array
        The ITRF coordinates of all of the baselines.
    bname : list
        The names of all of the baselines.
    """
    ant_itrf = get_itrf(
        latlon_center=(ct.OVRO_LAT*u.rad, ct.OVRO_LON*u.rad, ct.OVRO_ALT*u.m)
    ).loc[antennas]
    xx = np.array(ant_itrf['dx_m'])
    yy = np.array(ant_itrf['dy_m'])
    zz = np.array(ant_itrf['dz_m'])
    # Get uvw coordinates
    nants = len(antennas)
    nbls = (nants*(nants+1))//2
    blen = np.zeros((nbls, 3))
    bname = []
    k = 0
    for i in range(nants):
        for j in range(i, nants):
            blen[k, :] = np.array([
                xx[i]-xx[j],
                yy[i]-yy[j],
                zz[i]-zz[j]
            ])
            bname += ['{0}-{1}'.format(
                antennas[i],
                antennas[j]
            )]
            k += 1
    return blen, bname
