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


# set up defaults
candname = "260212slow"
cand_dir = "/dataz/dsa110/candidates/"+candname
msdir = "/media/ubuntu/ssd/localization_processing/"+candname
hdf5dir = cand_dir+"/Level2/calibration/"
#ra, dec = get_radec_T2_json(candname,cand_dir)
#print(ra, dec)
ra=344.65833*u.deg
dec=16.151389*u.deg
cal = generate_calibrator_source(candname, ra=ra, dec=dec)
duration = 15.*u.min # length of output
filelength = 5.*u.min # length of input files

# find files
def find_files(cname,cdir,mjd,refsb="sb00"):
    
    """
    cname (str): name of candidate.
    cdir (str): directory where candidate is archived (according to plan)
    refsb (str): sb used to find files (default sb00)
    Returns:
    list of files, date
    """

    # find files from candidate dir
    files = sorted(glob.glob(f"{cdir}/Level2/calibration/*_{refsb}.hdf5"))
    datetimes = [f.split("/")[-1][:19] for f in files]
    dates = np.unique([dt[:10] for dt in datetimes])
    transit_files = []    
    for dt in datetimes:
        if cal_in_datetime(dt,Time(mjd,format='mjd'), duration, filelength):            
            transit_files += [dt]

    return transit_files,dates[0]


# run main

mjd = 61083.89158

files,date = find_files(candname,cand_dir,mjd)
print("found files at date:",files,date)

# create ms
print("Creating ms...")
convert_calibrator_pass_to_ms(
    cal,
    date,
    files,
    msdir=msdir,
    hdf5dir=hdf5dir,
    refmjd=config.refmjd
)




    



    




