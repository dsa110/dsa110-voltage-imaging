"""Match sources between catalogs."""

import os
from typing import Tuple, Union

import numpy as np
import pandas
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.io import fits
import bdsf


class VLASSCat:
    """The VLASS catalog as a searchable dataframe."""

    def __init__(self, filepath: str = "/home/ubuntu/dana/vlass.csv"):
        self.filepath = filepath
        self._load()

    def _load(self):
        """Parse the vlass catalog."""
        vlass = pandas.read_csv(self.filepath)
        vlass.columns = [colname.strip() for colname in vlass.columns]
        vlass.drop(
            columns=[
                'Unnamed: 0', 'Isl_Total_flux', 'E_Isl_Total_flux', 'Isl_rms', 'Isl_mean',
                'Resid_Isl_rms', 'Resid_Isl_mean', 'e2_image', 'tile', 'img_coord',
                'e2_5x5_box_peak_mJy', 'e2_rms_mJy', 'e2_distance_to_img_edge_pix',
                'nlines', 'e2_flux_percentile', 'e1_image',
                'e1_5x5_box_peak_mJy', 'e1_rms_mJy', 'e1_distance_to_img_edge_pix',
                'e1_warning', 'corrected_e2_5x5_box_peak_mJy', 'systematic_fluxscale_corrections',
                'corrected_e2_rms_mJy', 'sys_residual_plus_err',
                'sys_residual_minus_err', 'true_E1_flux_mJy', 'E2_gtr_E1_quad_RMS',
                'p_e2_gtr_e1', 'true_corrected_E2_flux_mJy', 'p_e1_gtr_e2',
                'Modulation_index', 'variability_statistic', 'dist_5mJy', 'dist_10mJy',
                'dist_25mJy', 'dist_100mJy', 'dist_1Jy'],
            inplace=True)
        vlass = vlass[~((vlass["strict_stripe_warning"]) | (vlass["strict_sidelobe_warning"]))]
        vlass.drop(columns=["strict_stripe_warning", "strict_sidelobe_warning"], inplace=True)
        vlass['S_Code'] = [scode.strip() for scode in vlass['S_Code']]
        self._catalog = vlass

    def search(
            self, pointing: Union[Tuple[u.Quantity], SkyCoord], halfwidth: u.Quantity = 5*u.deg,
            code: str = "", min_flux: float = None) -> pandas.DataFrame:
        """Search for vlass sources.

        This currently does not use a search radius, but instead returns sources within a square
        region around `pointing`.

        Parameters
        ----------
        pointing :
            The coordinates around which to search.
        halfwidth :
            The maximum deviation in RA or Dec from `pointing` for which to return sources.
        code :
            If passed, restricts returned sources to one that match `code`. Can be "S", "M" or "C".
        min_flux :
            If passed, restricts returned sources to ones with Total_flux > `min_flux` mJy.
        """
        if not isinstance(pointing, SkyCoord):
            pointing = SkyCoord(*pointing)

        halfwidth_ra = halfwidth/np.cos(pointing.dec)
        ramin, ramax = [(pointing.ra + sign*halfwidth_ra).to_value(u.deg) for sign in [-1, 1]]
        decmin, decmax = [(pointing.dec + sign*halfwidth).to_value(u.deg) for sign in [-1, 1]]

        query_result = (
            (self._catalog['RA'] > ramin) & (self._catalog['RA'] < ramax) &
            (self._catalog['DEC'] > decmin) & (self._catalog['DEC'] < decmax))

        if code:
            query_result = query_result & (self._catalog['S_Code'] == code)
        if min_flux is not None:
            query_result = query_result & (self._catalog['Total_flux'] > min_flux)

        return self._catalog[query_result]

    @property
    def catalog(self) -> pandas.DataFrame:
        """Return a copy of the catalog."""
        return self._catalog.copy(deep=True)


class Image:
    """A radio image, stored in fits format."""

    def __init__(self, filepath: str):
        self._filepath = filepath
        self._load()
        self._sources = None

    def _load(self) -> None:
        """Read in an image, including the data and wcs coordinates."""
        with fits.open(self._filepath) as hdulist:
            hdu = hdulist[0]
            self._wcs = WCS(hdu.header)
            self._data = hdu.data[0, 0, ...]
        self._centre, self._f0, self._p0 = self._wcs.pixel_to_world(
            self._data.shape[0]//2, self._data.shape[1]//2, 0, 0)
        edge = self._wcs.pixel_to_world(0, 0, 0, 0)[0]
        self._radius = abs(self._centre.dec - edge.dec)
        self._datamean = np.mean(self._data)
        self._datastd = np.std(self._data)

    def find_sources(self, sparse=True, **kwargs) -> None:
        """Find sources in the image using pybdsf.

        If `sparse` is `True`, then sources close to bright sources are removed from the source
        list.
        kwargs are passed to pybdsf.
        """
        output = bdsf.process_image(self._filepath, **kwargs)
        output.write_catalog(format='csv', clobber=True)
        self.load_source_catalog(sparse)

    def load_source_catalog(self, sparse=True):
        """Load sources from a previously generated source catalog.

        If `sparse` is `True`, then sources close to bright sources are removed from the source
        list.
        """
        self._sources = read_bdsf_csv(f"{os.path.splitext(self._filepath)[0]}.pybdsm.gaul")
        if sparse:
            self.sparsify_source_list()

    def sparsify_source_list(self, minsep: u.Quantity = 10*u.arcminute):
        """Remove sources within `minsep` of a brighter source from the source catalog."""
        assert self._sources is not None
        sources = self._sources

        to_remove = []
        to_keep = []
        tol_dec = minsep.to_value(u.deg)
        tol_ra = tol_dec/np.cos(self._centre.dec)

        for i in range(len(sources)):
            if i not in to_remove + to_keep:
                ra = sources.iloc[i]["RA"]
                dec = sources.iloc[i]["DEC"]

                nearby_idxs = (
                    (sources["RA"] > (ra - tol_ra)) & (sources["RA"] < (ra + tol_ra)) &
                    (sources["DEC"] > (dec - tol_dec)) & (sources["DEC"] < (dec + tol_dec)))
                nearby = sources[nearby_idxs].sort_values(by='Total_flux', ascending=False)

                if len(nearby) == 0:
                    continue
                to_keep.extend(nearby.head(1).index)
                if len(nearby) > 1:
                    to_remove.extend(nearby.tail(len(nearby)-1).index)

        self._sources = sources.iloc[to_keep].drop_duplicates()
        self._sources.reset_index(inplace=True, drop=True)

    def show(self):
        """Plot the image, and return ax for overplotting additional elements."""
        fig, ax = plt.subplots(
            figsize=(16, 16),
            subplot_kw={'projection': dsaimage.wcs, 'slices': ('x', 'y', 0, 0)})
        im = ax.imshow(
            self._data, interpolation='nearest', cmap='gray', origin='lower',
            vmin=self._datamean - self._datastd * 2, vmax=self._datamean + self._datastd*10)
        return ax

    @property
    def filepath(self) -> str:
        """The path to the fits image on disk."""
        return self._filepath

    @property
    def data(self) -> np.ndarray:
        """The image as a numpy array."""
        return self._data

    @property
    def wcs(self) -> WCS:
        """The WCS coordinates of the image."""
        return self._wcs

    @property
    def centre(self) -> SkyCoord:
        """The coordinates at the centre of the image."""
        return self._centre

    @property
    def f0(self) -> u.Quantity:
        """The central frequency of the image."""
        return self._f0

    @property
    def p0(self) -> u.Quantity:
        """The central polarization channel of the image."""
        return self._p0

    @property
    def radius(self) -> u.Quantity:
        """The radius of the image."""
        return self._radius

    @property
    def sources(self) -> pandas.DataFrame:
        """The sources in the image."""
        if self._sources is None:
            try:
                self.load_source_catalog()
                print("Loaded source catalog from disk")
            except FileNotFoundError:
                self.find_sources()
                print("Generated new source catalog")

        return self._sources


def read_bdsf_csv(filepath: str) -> pandas.DataFrame:
    """Read in sources from a csv generated by bdsf."""
    sources = pandas.read_csv(filepath, header=4)
    sources.columns = [colname.strip() for colname in sources.columns]
    sources['S_Code'] = [scode.strip() for scode in sources['S_Code']]
    return sources


def match_catalogs(dsa: pandas.DataFrame, reference: pandas.DataFrame) -> pandas.DataFrame:
    """Match sources from two catalogs.

    Beware: the ids returned are the ilocs, not indexs."""
    Match = namedtuple('Match', 'dsa_id ref_id sep_2d_arcsecond ra_offset_deg dec_offset_deg ref_code')
    dsa_coords = SkyCoord(dsa['RA']*u.deg, dsa['DEC']*u.deg)
    ref_coords = SkyCoord(reference['RA']*u.deg, reference['DEC']*u.deg)
    matches = []
    for dsa_id, dsa_coord in enumerate(dsa_coords):
        ref_id, sep2d, _ = match_coordinates_sky(dsa_coord, ref_coords)

        delta_ra = (ref_coords[ref_id].ra - dsa_coord.ra).to_value(u.arcsec)/np.cos(dsa_coord.dec)
        delta_dec = (ref_coord[ref_id].dec - dsa_coord.dec).to_value(u.arcsec)

        matches.append(Match(
            dsa_id, ref_id, sep2d[0].to_value(u.arcsecond), delta_ra, delta_dec,
            reference.iloc[ref_id]['S_Code']))

    return pandas.DataFrame(matches, columns=Match._fields)
