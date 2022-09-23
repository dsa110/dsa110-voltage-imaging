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

# set up defaults
#candname = "220411aabk"
candname = sys.argv[1]
cand_dir = "/dataz/dsa110/candidates/"+candname
msdir = cand_dir+"/Level2/calibration/"
hdf5dir = cand_dir+"/Level2/calibration/"
cal = generate_calibrator_source(candname, ra=None, dec=None)
duration = 5.*u.min # length of output
filelength = 5.*u.min # length of input files

# find files
def find_files_with_T2_json(cname,cdir,refsb="sb00"):
    
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
    mjd = canddata['mjds']
    files = sorted(glob.glob(f"{cdir}/Level2/calibration/*_{refsb}.hdf5"))
    datetimes = [f.split("/")[-1][:19] for f in files]
    dates = np.unique([dt[:10] for dt in datetimes])
    transit_files = []
    for dt in datetimes:
        if cal_in_datetime(dt,Time(mjd,format='mjd'), duration, filelength):            
            transit_files += [dt]

    return transit_files,dates[0]


# run main

files,date = find_files_with_T2_json(candname,cand_dir)
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
    
    



    




