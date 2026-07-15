import numpy as np
from casatasks import *
import sys
import glob
import os
import casatools as cc
from os import path
from pathlib import Path

candname = sys.argv[1]
flags = sys.argv[2]
fflags = flags.split(',')
flagstr = f"0:{fflags[0]}"
for fl in fflags[1:]:
    flagstr = f"{flagstr},0:{fl}"

wdir = f"/dataz/dsa110/localization_processing/{candname}/"

# calibrate voltage ms
frbms = f"{wdir}/{candname}.ms"
oms = f"{wdir}/voltage.ms"
os.system(f"rm -rf {oms}*")
split(vis=frbms,outputvis=oms,datacolumn="data")

flagdata(vis=oms,spw=flagstr)
flagdata(vis=oms,antenna="10")

os.system(f"mkdir -p {wdir}/icrs_image")
os.system(f"mkdir -p {wdir}/image")

my_file = Path(f"{wdir}/rfc_idict_more.npz")

if not my_file.is_file():
    sys.exit(1)
    
idict = np.load(f"{wdir}/rfc_idict_more.npz",allow_pickle=True)["idict"].item()
for src in idict['jname']:
    try:
        if os.path.exists(f"{wdir}/icrs_image/{src}-image.fits"):
            print(f"Image already made for {src}")
        else:
            myms = glob.glob(f"{wdir}/*{src}.ms")[0]
            print(f"Calibrating {myms}")
            applycal(vis=myms,gaintable=[f"{wdir}/cal.B0"])
            print(f"Flagging {myms}")
            flagdata(vis=myms,spw=flagstr)
            flagdata(vis=myms,antenna="10")
            split(vis=myms,outputvis=f"{wdir}/icrs_image/tmp.ms",width=64)
            os.system(f"cd {wdir}/icrs_image; wsclean -weight briggs 0.5 -name {src} -size 1200 1200 -scale 0.2asec -minuv-l 300 -niter 20 -auto-threshold 4.0 -theoretic-beam -nwlayers 1 tmp.ms; cd -")
            os.system(f"rm -rf {wdir}/icrs_image/tmp.ms")
    except:
        print(f"Could not find {src}")

