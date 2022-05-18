"""Match sources between catalogs."""

import os
from typing import Tuple, Union, List

import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Arrow
import numpy as np
import pandas
import astropy.units as u
from astropy.coordinates import SkyCoord, match_coordinates_sky
from astropy.wcs import WCS
from astropy.io import fits
import bdsf


__all__ = [
    "VLASSCat", "Image", "SourceCatalog", "match_catalogs", "plot_catalogs",
    "plot_offset_direction", "plot_matched_sources", "plot_offsets"]


class VLASSCat:
    """The VLASS catalog as a searchable dataframe."""

    def __init__(self, filepath: Union[str, "Path"] = "/home/ubuntu/dana/vlass.csv"):
        self.filepath = str(filepath)
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


class SourceCatalog:
    """Sources found in an Image."""

    def __init__(self, image: "Image", catalogpath: Union[str, "Path"] = ""):
        self._imagepath = image.filepath
        if catalogpath:
            self._catalogpath = str(catalogpath)
        else:
            self._catalogpath = f"{os.path.splitext(self._imagepath)[0]}.pybdsm.gaul"
        self._sources = None

    def generate(self, **kwargs):
        """Find sources in the image using pybdsf.

        If `sparse` is `True`, then sources close to bright sources are removed from the source
        list.
        kwargs are passed to pybdsf.
        """
        output = bdsf.process_image(self._imagepath, **kwargs)
        output.write_catalog(format='csv', outfile=self._catalogpath, clobber=True)
        if not os.path.exists(self._catalogpath):
            self._write_empty_catalog()

    def load(self, sparse=True):
        """Load sources from a previously generated source catalog.

        If `sparse` is `True`, then sources close to bright sources are removed from the source
        list.
        """

        if (not os.path.exists(self._catalogpath) or
                os.path.getmtime(self._catalogpath) < os.path.getmtime(self._imagepath)):
            self.generate()

        self._sources = read_bdsf_csv(self._catalogpath)
        if sparse:
            self.sparsify()

    def sparsify(self, minsep: u.Quantity = 10*u.arcminute):
        """Remove sources within `minsep` of a brighter source from the source catalog."""
        # TODO: Query sources within a circle instead of a square
        assert self._sources is not None
        sources = self._sources
        centre_dec = np.mean(sources['DEC'])

        to_remove = []
        to_keep = []
        tol_dec = minsep.to_value(u.deg)
        tol_ra = tol_dec/(np.cos(centre_dec*u.deg)).to_value(u.dimensionless_unscaled)

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

    def _write_empty_catalog(self):
        with open(self._catalogpath, 'w', encoding='utf-8') as f:
            f.write((
                f"# Gaussian list for {self._imagepath}\n"
                "# Generated by PyBDSM version 1.10.1\n"
                f"# Reference frequency of the detection (\"ch0\") image: Unknown\n"
                "# Equinox : 2000.0\n"
                "\n"
                "# Gaus_id, Isl_id, Source_id, Wave_id, RA, E_RA, DEC, E_DEC, Total_flux, "
                "E_Total_flux, Peak_flux, E_Peak_flux, Xposn, E_Xposn, Yposn, E_Yposn, Maj, "
                "E_Maj, Min, E_Min, PA, E_PA, Maj_img_plane, E_Maj_img_plane, Min_img_plane, "
                "E_Min_img_plane, PA_img_plane, E_PA_img_plane, DC_Maj, E_DC_Maj, DC_Min, "
                "E_DC_Min, DC_PA, E_DC_PA, DC_Maj_img_plane, E_DC_Maj_img_plane, "
                "DC_Min_img_plane, E_DC_Min_img_plane, DC_PA_img_plane, E_DC_PA_img_plane, "
                "Isl_Total_flux, E_Isl_Total_flux, Isl_rms, Isl_mean, Resid_Isl_rms, "
                "Resid_Isl_mean, S_Code\n"))

    @property
    def sources(self) -> pandas.DataFrame:
        """The sources in the image."""
        if self._sources is None:
            self.load()
        return self._sources

    @property
    def catalogpath(self) -> str:
        """The path to the catalog containing found sources."""
        return self._catalogpath

    @catalogpath.setter
    def catalogpath(self, value: Union[str, "Path"]) -> None:
        """The path to the catlog containing found sources."""
        self._catalogpath = str(value)
        self._sources = None


class Image:
    """A radio image, stored in fits format."""

    def __init__(self, filepath: Union[str, "Path"]):
        self._filepath = str(filepath)
        self.catalog = SourceCatalog(self)
        self._load()

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

    def show(self):
        """Plot the image, and return ax for overplotting additional elements."""
        fig, ax = plt.subplots(
            figsize=(16, 16),
            subplot_kw={'projection': self._wcs, 'slices': ('x', 'y', 0, 0)})
        _im = ax.imshow(
            self._data, interpolation='nearest', cmap='gray', origin='lower',
            vmin=self._datamean - self._datastd * 2, vmax=self._datamean + self._datastd*10)
        return fig, ax

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
        return self.catalog.sources

def read_bdsf_csv(filepath: Union[str, "Path"]) -> pandas.DataFrame:
    """Read in sources from a csv generated by bdsf."""
    sources = pandas.read_csv(str(filepath), header=4)
    sources.columns = [colname.strip() for colname in sources.columns]
    sources['S_Code'] = [scode.strip() for scode in sources['S_Code']]
    return sources


def match_catalogs(dsa: pandas.DataFrame, reference: pandas.DataFrame) -> pandas.DataFrame:
    """Match sources from two catalogs.

    Beware: the ids returned are the ilocs, not indexs."""
    dsa_coords = SkyCoord(dsa['RA']*u.deg, dsa['DEC']*u.deg)
    ref_coords = SkyCoord(reference['RA']*u.deg, reference['DEC']*u.deg)
    ref_ids, sep2ds, _ = match_coordinates_sky(dsa_coords, ref_coords)

    sep2ds = [sep2d.to_value(u.arcsec) for sep2d in sep2ds]
    delta_ras = (
        (ref_coords[ref_ids].ra - dsa_coords.ra).to_value(u.arcsec) *
        np.cos(dsa_coords.dec.to_value(u.rad)))
    delta_decs = (ref_coords[ref_ids].dec - dsa_coords.dec).to_value(u.arcsec)
    dsa_ids = range(len(dsa_coords))
    ref_codes = reference.iloc[ref_ids]['S_Code']

    return pandas.DataFrame(
        list(zip(dsa_ids, ref_ids, sep2ds, delta_ras, delta_decs, ref_codes)),
        columns=['dsa_id', 'ref_id', 'sep2d_asec', 'ra_offset_asec', 'dec_offset_asec', 'ref_code'])

def plot_catalogs(
        dsaimage: Image, catalogs: List[Tuple[SkyCoord, int, str]]
) -> Tuple["matplotlib.fig", "matplotlib.axes.Axes"]:
    """Plot sources from two catalogs.

    `catalogs` is a list of tuples, which are
    (catalog coordinates, size of circle to draw, color of circle to draw)
    """
    fig, ax = dsaimage.show()
    for catalog, size, colour in catalogs:
        for i, coord in enumerate(catalog):
            x, y, *_ = dsaimage.wcs.world_to_pixel(coord, dsaimage.f0, dsaimage.p0)
            if 0 <= x < dsaimage.data.shape[0] and 0 <= y < dsaimage.data.shape[1]:
                e = Circle(xy=(x, y), radius=size)
                e.set_facecolor('none')
                e.set_edgecolor(colour)
                ax.add_artist(e)
                ax.annotate(i, (x, y), color=colour)
    return fig, ax

def plot_offset_direction(
        dsaimage: Image, coords: SkyCoord, ra_offsets: List[float], dec_offsets: List[float]
) -> Tuple["matplotlib.fig", "matplotlib.axes.Axes"]:
    """Plot measured offsets on an image."""
    fig, ax = dsaimage.show()
    for (coord, ra_offset, dec_offset) in zip(coords, ra_offsets, dec_offsets):
        x, y, *_ = dsaimage.wcs.world_to_pixel(coord, dsaimage.f0, dsaimage.p0)
        if 0 <= x < dsaimage.data.shape[0] and  0 <= y < dsaimage.data.shape[1]:
            e = Arrow(x=x, y=y, dx=ra_offset*200, dy=dec_offset*200, width=200)
            e.set_facecolor('white')
            e.set_edgecolor('white')
            ax.add_artist(e)
    return fig, ax

def plot_matched_sources(
        dsaimage: Image, catalogs: List[Tuple[SkyCoord, int, str]], matched: List[int],
        matched_size: int = 75, matched_colour: str = "yellow"
) -> Tuple["matplotlib.fig", "matplotlib.axes.Axes"]:
    """Plot two catalogs, as well as matches between them.

    `matched` are indices in the first catalog for which matches were found.
    `catalogs` is a list of tuples, which are
    (catalog coordinates, size of circle to draw, colour of circle to draw)
    """
    fig, ax = plot_catalogs(dsaimage, catalogs)
    coords = catalogs[0][0]
    for i in matched:
        x, y, *_ = dsaimage.wcs.world_to_pixel(coords[i], dsaimage.f0, dsaimage.p0)
        if 0 <= x < dsaimage.data.shape[0] and 0 <= y < dsaimage.data.shape[1]:
            e = Circle(xy=(x, y), radius=matched_size)
            e.set_facecolor('none')
            e.set_edgecolor(matched_colour)
            ax.add_artist(e)
    return fig, ax

def plot_offsets(
        ra_offsets: List[float], dec_offsets: List[float], limit: int = 4
) -> Tuple["matplotlib.fig", "matplotlib.axes.Axes"]:
    """Create a scatter plots of the offsets."""
    fig, ax = plt.subplots(1, 1, figsize=(6, 6))
    ax.scatter(ra_offsets, dec_offsets)
    ax.set_xlim(-1*limit, limit)
    ax.set_ylim(-1*limit, limit)
    ax.axhline(0, color='grey')
    ax.axvline(0, color='Grey')
    ax.set_xlabel('RA offset (arcsec)')
    ax.set_ylabel('Dec offset (arcsec)')
    return fig, ax
