"""Test the update of a template ms."""

from shutil import copytree
from typing import Union
from dsaT3.update_template import *

SINGLE_CORR_TEMPLATE = '/media/ubuntu/data/dsa110/imaging/template_corr08.ms'
SINGLE_CORR_MS = '/media/ubuntu/data/dsa110/imaging/data_corr08.ms'
SINGLE_CORR_UVH5_SAME = '/media/ubuntu/ssd/data/data_sameobs_corr08.hdf5'
SINGLE_CORR_UVH5_DIFF = '/media/ubuntu/ssd/data/data_diffobs_corr08.hdf5'

TOL = 1e-5

def test_update_visibilities_same(tmpdir: str, tol: float=TOL) -> None:
    """Test that visibilities are updated succesfully."""
    template_filepath = f'{tmpdir}/template_testvissame.ms'
    uvh5_filepaths = [SINGLE_CORR_UVH5_SAME]

    copytree(SINGLE_CORR_TEMPLATE, template_filepath)
    update_visibilities(template_filepath, uvh5_filepaths, n_corr_nodes=1)

    UV = UVData()
    UV.read(uvh5_filepaths[0], file_type='uvh5')
    
    uvh5_vis = UV.data_array.squeeze(1)
    uvh5_flag = UV.flag_array.squeeze(1)
    if np.median(np.diff(UV.freq_array.squeeze(0))) < 0:
        uvh5_vis = uvh5_vis[:, ::-1, :]
        uvh5_flag = uvh5_flag[:, ::-1, :]

    with table(template_filepath) as tb:
        assert arrays_equal(np.array(tb.DATA[:]), uvh5_vis, tol)
        assert np.all(np.array(tb.FLAG[:]) == uvh5_flag)

        with table(SINGLE_CORR_MS) as tb2:
            assert arrays_equal(np.array(tb.DATA[:]), np.array(tb2.DATA[:]), tol)
            assert np.all(np.array(tb.FLAG[:]) == np.array(tb2.FLAG[:]))

def test_update_visibilities_diff(tmpdir: str, tol: float=TOL) -> None:
    """Test that visibilities are updated succesfully."""
    template_filepath = f'{tmpdir}/template_testvisdiff.ms'
    uvh5_filepaths = [SINGLE_CORR_UVH5_DIFF]

    copytree(SINGLE_CORR_TEMPLATE, template_filepath)
    update_visibilities(template_filepath, uvh5_filepaths, n_corr_nodes=1)

    UV = UVData()
    UV.read(uvh5_filepaths[0], file_type='uvh5')
    
    uvh5_vis = UV.data_array.squeeze(1)
    uvh5_flag = UV.flag_array.squeeze(1)
    if np.median(np.diff(UV.freq_array.squeeze(0))) < 0:
        uvh5_vis = uvh5_vis[:, ::-1, :]
        uvh5_flag = uvh5_flag[:, ::-1, :]

    with table(template_filepath) as tb:
        assert arrays_equal(np.array(tb.DATA[:]), uvh5_vis, tol)
        assert np.all(np.array(tb.FLAG[:]) == uvh5_flag)

def test_update_metadata_same(tmpdir: str, tol: float=TOL) -> None:
    """Test that metadata is updated successfully."""
    # TODO: use a data path from a different observation
    template_filepath = f'{tmpdir}/template_testmdsame.ms'
    uvh5_filepath = SINGLE_CORR_UVH5_SAME

    copytree(SINGLE_CORR_TEMPLATE, template_filepath)
    update_metadata(template_filepath, uvh5_filepath)

    UV = UVData()
    UV.read(uvh5_filepath, file_type='uvh5')

    assert_times_match(UV, template_filepath, tol)
    assert_uvws_match(UV, template_filepath, tol)
    assert_directions_match(UV, template_filepath, tol)

def test_update_metadata_diff_obs(tmpdir: str, tol: float=TOL) -> None:
    """Test that metadata is updated successfully using a different 
    observation than was used to generate the template.
    """
    # TODO: use a data path from a different observation
    template_filepath = f'{tmpdir}/template_testmddiff.ms'
    uvh5_filepath = SINGLE_CORR_UVH5_DIFF

    copytree(SINGLE_CORR_TEMPLATE, template_filepath)
    update_metadata(template_filepath, uvh5_filepath)

    UV = UVData()
    UV.read(uvh5_filepath, file_type='uvh5')

    assert_times_match(UV, template_filepath, tol)
    assert_uvws_match(UV, template_filepath, tol)
    assert_directions_match(UV, template_filepath, tol)

def assert_times_match(UV, template_filepath, tol) -> None:
    """Assert that times betwen the template file and the UV object match."""
    time_offset = -0.00011921
    tobs_mjds = convert_jd_to_mjds(UV.time_array)
    tstart = tobs_mjds[0]+time_offset
    tstart_source = tobs_mjds[0]

    with table(template_filepath) as tb:
        assert arrays_equal(np.array(tb.TIME[:]), tobs_mjds, tol)
        assert arrays_equal(np.array(tb.TIME_CENTROID[:]), tobs_mjds, tol)

    with table(f'{template_filepath}/FEED') as tb:
        assert arrays_equal(np.array(tb.TIME[:]), tstart, tol)

    with table(f'{template_filepath}/FIELD') as tb:
        assert arrays_equal(np.array(tb.TIME[:]), tstart, tol)

    with table(f'{template_filepath}/OBSERVATION') as tb:
        assert arrays_equal(np.array(tb.TIME_RANGE[:]), tstart, tol)

    with table(f'{template_filepath}/SOURCE') as tb:
        assert arrays_equal(np.array(tb.TIME[:]), tstart_source, tol)

def assert_uvws_match(UV, template_filepath, tol) -> None:
    """Assert that uvw's between the template file and the UV object match."""
    with table(template_filepath) as tb:
        assert arrays_equal(np.array(tb.UVW[:]), UV.uvw_array, tol)

def assert_directions_match(UV, template_filepath, tol) -> None:
    """Update the directions in the template ms with the true direction."""
    ra_rad, dec_rad = get_pointing(UV)

    with table(f'{template_filepath}/FIELD') as tb:
        for i, value in enumerate([ra_rad, dec_rad]):
            print(value)
            assert arrays_equal(np.array(tb.DELAY_DIR[:])[..., i], value, tol)
            assert arrays_equal(np.array(tb.PHASE_DIR[:])[..., i], value, tol)
            assert arrays_equal(np.array(tb.REFERENCE_DIR[:])[..., i], value, tol)

    with table(f'{template_filepath}/SOURCE') as tb:
        for i, value in enumerate([ra_rad, dec_rad]):
            assert arrays_equal(np.array(tb.DIRECTION[:])[..., i], value, tol)

def arrays_equal(array: np.array, value: Union[np.array, float], tol):
    """Check that two arrays or an array and a value are equal.

    If `value` is an array, then `array` and `value` must have the same shape.
    """
    if isinstance(value, np.ndarray):
        if not array.shape == value.shape:
            return False
    if not np.all(np.abs(array-value) < tol):
        return False

    return True
