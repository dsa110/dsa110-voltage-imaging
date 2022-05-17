import astropy.units as u
import pandas
from typing import Tuple
from astropy.coordinates import SkyCoord

class VLASSCat:
    
    def __init__(self, filepath: str = "/home/ubuntu/dana/vlass.csv"):
        self._load(filepath)

    def _load(self, filepath: str):
        """Parse the vlass catalog."""
        vlass = pandas.read_csv(filepath)
        vlass.columns = [colname.strip() for colname in vlass.columns]
        vlass.drop(
            columns=[
                'Unnamed: 0', 'Isl_Total_flux', 'E_Isl_Total_flux', 'Isl_rms', 'Isl_mean', 'Resid_Isl_rms',
                'Resid_Isl_mean', 'e2_image', 'tile', 'img_coord',
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

    def search(self, pointing: Tuple["Quantity"], halfwidth: "Quantity" = 5*u.deg, code: str = ""):
        """Search for compact vlass sources."""
        if not isinstance(pointing, SkyCoord):
            pointing = SkyCoord(*pointing)
        
        ramin, ramax = [(pointing.ra + sign*halfwidth).to_value(u.deg) for sign in [-1, 1]]
        decmin, decmax = [(pointing.dec + sign*halfwidth).to_value(u.deg) for sign in [-1, 1]]

        query_result = (
            (self._catalog['RA'] > ramin) & (self._catalog['RA'] < ramax) & 
            (self._catalog['DEC'] > decmin) & (self._catalog['DEC'] < decmax))

        if code:
            query_result = query_result & (self._catalog['S_Code'] == code)
        
        return self._catalog[query_result]

    @property
    def catalog(self):
        """Return a copy of the catalog."""
        return self._catalog.copy(deep=True)
    

# class SourceFinder:
#     def __init__(self):
#         pass


# class CatMatcher:
#     def __init__(self):
#         pass

def read_bdsf_csv(filepath: str) -> pandas.DataFrame:
    sources = pandas.read_csv(filepath, header=4)
    sources.columns = [colname.strip() for colname in sources.columns]
    sources['S_Code'] = [scode.strip() for scode in sources['S_Code']]
    return sources
