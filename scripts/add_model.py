from dsacalib.ms_io import convert_calibrator_pass_to_ms
from dsacalib.utils import generate_calibrator_source
from dsacalib.routines import cal_in_datetime
import astropy.units as u
import numpy as np
import json
import glob
from astropy.time import Time
import dsacalib.config as configuration
config = configuration.Configuration()
import sys
from astropy.coordinates import SkyCoord
from dsacalib.preprocess import read_nvss_catalog
import pandas as pd
pd.options.mode.chained_assignment = None
df = read_nvss_catalog()
import casatools as cc
import casatasks as cst
import os
from dsautils import coordinates
import sys

# function to add model to ms, returning model components
def add_model_to_ms(msdata,flux_limit=50.,radius_limit=2.,complist_name='mycomplist.cl',ra_pt=100.0,dec_pt=71.66):
    """
    Adds NVSS components to model_data column (overwrites). 
    Inputs:
       msdata (str): path to measurement set [no default]
       flux_limit (float, mJy): minimum flux density of model components [50 mJy]
       radius_limit (float, deg): max radius of components from phase center [2 deg]
       complist_name (str): path/name of component list file [mycomplist.cl]
       ra_pt (float, deg) [default 100.0]
       dec_pt (float, deg) [default 71.66]
    Returns:
        ra, dec, flux of components. 
    """    

    # do cross match
    ra = df.loc[:,'ra']; dec = df.loc[:,'dec']; c = SkyCoord(ra,dec,frame='icrs',unit='deg')
    pos = SkyCoord([ra_pt],[dec_pt],frame='icrs',unit='deg')
    print(pos)
    idx, d2d, d3d = c.match_to_catalog_sky(pos)
    
    # select sources from NVSS
    flux = df.loc[:,'flux_20_cm']
    mat = (df.iloc[np.where((d2d.deg<radius_limit) & (flux>flux_limit))])
    
    # extract columns, accounting for extended sources
    m_ra = mat.loc[:,'ra']; m_h = np.floor(m_ra/15.).astype('int'); m_m = np.floor(60.*(m_ra/15.-np.floor(m_ra/15.))).astype('int'); m_s = 60.*(60.*(m_ra/15.-np.floor(m_ra/15.))-1.*m_m )
    m_dec = mat.loc[:,'dec']; m_dd = np.floor(m_dec).astype('int'); m_dm = np.floor(60.*(m_dec-np.floor(m_dec))).astype('int'); m_ds = 60.*(60.*(m_dec-np.floor(m_dec))-1.*m_dm )
    m_flux = mat.loc[:,'flux_20_cm']
    m_maj = mat.loc[:,'major_axis']
    m_min = mat.loc[:,'minor_axis']
    m_pa = mat.loc[:,'position_angle']*np.pi/180.

    m_maj[np.isnan(m_pa)] = 0.
    m_min[np.isnan(m_pa)] = 0.
    m_min[m_maj<30.] = 0.
    m_pa[m_maj<30.] = 0.
    m_maj[m_maj<30.] = 0.
    m_pa[np.isnan(m_pa)] = 0.
    
    # make component list
    os.system(f"if [ -e {complist_name} ]; then rm -rf {complist_name}; fi")
    a = cc.componentlist()
    for i in np.arange(len(m_ra)):
        direct = "J2000 %2dh%dm%.2fs +%dd%dm%.2fs"%(m_h[i],m_m[i],m_s[i],m_dd[i],m_dm[i],m_ds[i])
        print(direct)
        if m_maj[i]==0.:
            a.addcomponent(shape="Point",dir=direct,flux=m_flux[i]*0.001,fluxunit='Jy',freq='1.405GHz')
        else:
            majj = "%.1farcsec"%(m_maj[i])
            minn = "%.1farcsec"%(m_min[i])
            posang = "%.1fdeg"%(m_pa[i])
            a.addcomponent(shape="Gaussian",dir=direct,flux=m_flux[i]*0.001,fluxunit='Jy',freq='1.405GHz',majoraxis=majj,minoraxis=minn,positionangle=posang)
    a.rename(complist_name)
    a.close()
    print(f"Made component list {complist_name}")
    
    # ft into model column
    cst.ft(vis=msdata,complist=complist_name,spw='0',usescratch=True)
    
    return m_ra,m_dec,m_flux

# get RA and DEC to phase to
def get_radec_T2_json(cname,cdir):
    """
    cname (str): name of candidate.
    cdir (str): directory where candidate is archived (according to plan)
    Returns:
    ra, dec
    """

    # load in T2 json
    try:
        canddata = json.load(open(cdir+"/Level2/voltages/T2_"+cname+".json"))
    except:
        print("Couldnt open T2 json.")
        return None

    # get ra and dec
    canddata = canddata[cname]
    mjd = canddata['mjds']
    ra, dec = coordinates.get_pointing(ibeam=128,obstime=Time(mjd,format='mjd'))
    return ra, dec
    
ra=282.44
dec=71.66


print("Adding model...")
ms = sys.argv[1]
cst.delmod(ms)
print(f"Ready with {ms}")
my_ra,my_dec,my_flux = add_model_to_ms(ms,ra_pt=ra,dec_pt=dec)




    



    




