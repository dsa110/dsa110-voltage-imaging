from dsacalib.ms_io import convert_calibrator_pass_to_ms
from dsacalib.utils import generate_calibrator_source
from dsacalib.routines import cal_in_datetime
from dsavim.utils import rfc_lookup
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

# get frbpos (M8 archive; replaces the T3 T2_<cand>.json reader)
def frbpos_from_archive(cname):
    from dsavim.metadata import load_burst
    burst = load_burst(cname)
    return SkyCoord(burst.ra_deg*u.deg, burst.dec_deg*u.deg), burst.t_peak.mjd

# find files
def find_transit_files(cname,cdir,dt=0.,refsb="sb00"):
    
    """
    cname (str): name of candidate.
    cdir (str): directory where candidate is archived (according to plan)
    refsb (str): sb used to find files (default sb00)
    Returns:
    list of files, date
    """
    
    mjd = cand_mjd + dt
    files = sorted(glob.glob(f"{cdir}/calibration/*_{refsb}.hdf5"))
    datetimes = [f.split("/")[-1][:19] for f in files]
    dates = np.unique([dt[:10] for dt in datetimes])
    transit_files = []
    for dt in datetimes:
        if cal_in_datetime(dt,Time(mjd,format='mjd'), duration, filelength):            
            transit_files += [dt]

    return transit_files,dates[0]


# set up defaults
candname = sys.argv[1]
cand_dir = "/dataz/dsa110/candidates/"+candname
hdf5dir = cand_dir+"/calibration/"
duration = 5.*u.min # length of output
filelength = 5.*u.min # length of input files

# generate cal source
frbpos, cand_mjd = frbpos_from_archive(candname)
idict = rfc_lookup(frbpos,dra=8.5,ddec=1.3,thresh=0.02,findNearest=False,findBrightest=True,nBrightest=10)
print(idict)

# make directory and save idict
odir = "/media/ubuntu/ssd/localization_processing/"+candname+"/"
msdir = odir
np.savez(odir+"rfc_idict.npz",idict=idict)


for i in np.arange(len(idict['sep'])):
    cal = generate_calibrator_source(idict['jname'][i],ra=idict['position'][i].ra,dec=idict['position'][i].dec)
    dt = ((idict['position'][i].ra.deg - frbpos.ra.deg)/360.) # in days
    files,date = find_transit_files(candname,cand_dir,dt=dt)
    print("found files at date:",files,date)

    # create ms
    print(f"Creating ms {i+1} of {len(idict['sep'])}...")
    convert_calibrator_pass_to_ms(
        cal,
        date,
        files,
        msdir=msdir,
        hdf5dir=hdf5dir,
        refmjd=config.refmjd
    )


    



    




