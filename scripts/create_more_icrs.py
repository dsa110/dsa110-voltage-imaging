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

# get frbpos
def frbpos_from_json(cname,cdir):

    # load in T2 json
    try:
        canddata = json.load(open(cdir+"/Level2/voltages/T2_"+cname+".json"))
    except:
        print("Couldnt open T2 json.")
        return None

    # find files from candidate dir
    canddata = canddata[cname]
    ra = canddata["ra"]
    dec = canddata["dec"]

    return SkyCoord(ra*u.deg,dec*u.deg)

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
inspect = False
if len(sys.argv)>2:
    inspect=True
cand_dir = "/dataz/dsa110/candidates/"+candname
hdf5dir = cand_dir+"/Level2/calibration/"
duration = 5.*u.min # length of output
filelength = 5.*u.min # length of input files
odir = "/media/ubuntu/data/localization_processing/"+candname+"/"
msdir = odir


# generate cal sources
frbpos = frbpos_from_json(candname,cand_dir)
idict = rfc_lookup(frbpos,dra=8.5,ddec=1.5,thresh=0.01,findNearest=False,findBrightest=True,nBrightest=24)
#print(idict)

# find existing sources
mss = glob.glob(f"{odir}/*_J*.ms")
existing_ms = []
for fl in mss:
    existing_ms.append(fl[-13:-3])
made_file = []
for src in idict['jname']:
    if src in existing_ms:
        made_file.append(True)
    else:
        made_file.append(False)

myn = 13
nfiles = 0
make_file = []
for i in np.arange(len(made_file)):

    dt = ((idict['position'][i].ra.deg - frbpos.ra.deg)/360.) # in days
    files,date = find_files_with_T2_json(candname,cand_dir,dt=dt)
    if len(files)>0:
        print(i,'.. CAN DO:',idict['jname'][i],idict['flux'][i],idict['sep'][i],made_file[i])
        if made_file[i] is True:
            nfiles += 1
            make_file.append(False)
        if made_file[i] is False:
            if nfiles<myn:
                print('Will process above...')
                make_file.append(True)
                nfiles += 1
            else:
                make_file.append(False)
    else:
        print(i,'.. NOT IN DATA:',idict['jname'][i],idict['flux'][i],idict['sep'][i],made_file[i])
        make_file.append(False)

if inspect:
    sys.exit(1)
    
np.savez(odir+"rfc_idict_more.npz",idict=idict)


for i in np.arange(len(idict['sep'])):
    cal = generate_calibrator_source(idict['jname'][i],ra=idict['position'][i].ra,dec=idict['position'][i].dec)
    dt = ((idict['position'][i].ra.deg - frbpos.ra.deg)/360.) # in days
    files,date = find_files_with_T2_json(candname,cand_dir,dt=dt)
    print("found files at date:",files,date)

    # create ms

    if make_file[i] is True:
        print(f"Creating ms...")
        convert_calibrator_pass_to_ms(
            cal,
            date,
            files,
            msdir=msdir,
            hdf5dir=hdf5dir,
            refmjd=config.refmjd
        )
        


    



    




