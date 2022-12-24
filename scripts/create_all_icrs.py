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

# find files
def find_files_with_T2_json(cname,cdir,dt=0.,refsb="sb00"):
    
    """
    cname (str): name of candidate.
    cdir (str): directory where candidate is archived (according to plan)
    refsb (str): sb used to find files (default sb00)
    Returns:
    list of files, date
    """
    
    # load in T2 json
    try:
        canddata = json.load(open(cdir+"/Level2/voltages/T2_"+cname+".json"))
    except:
        print("Couldnt open T2 json.")
        return None

    # find files from candidate dir
    canddata = canddata[cname]
    mjd = canddata['mjds'] + dt
    files = sorted(glob.glob(f"{cdir}/Level2/calibration/*_{refsb}.hdf5"))
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
hdf5dir = cand_dir+"/Level2/calibration/"
duration = 5.*u.min # length of output
filelength = 5.*u.min # length of input files

# generate cal source
ra = sys.argv[2]
dec = sys.argv[3]
frbpos = SkyCoord(ra,dec,unit=(u.hourangle, u.deg))
idict = rfc_lookup(frbpos,dra=8.5,ddec=1.0,thresh=0.04,findNearest=False)

# make directory and save idict
odir = "/media/ubuntu/data/vikram/localization_test/"+candname+"/"
msdir = odir
os.system("mkdir -p "+odir)
np.savez(odir+"idict.npz",idict=idict)

print(idict)

for i in np.arange(len(idict['sep'])):
    cal = generate_calibrator_source(idict['jname'][i],ra=idict['position'][i].ra,dec=idict['position'][i].dec)
    dt = ((idict['position'][i].ra.deg - frbpos.ra.deg)/360.) # in days
    files,date = find_files_with_T2_json(candname,cand_dir,dt=dt)
    print("found files at date:",files,date)

    # create ms
    print(f"Creating ms {i+1} of 2...")
    convert_calibrator_pass_to_ms(
        cal,
        date,
        files,
        msdir=msdir,
        hdf5dir=hdf5dir,
        refmjd=config.refmjd
    )


    



    




