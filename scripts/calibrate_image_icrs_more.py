import numpy as np
from casatasks import *
import sys
import glob
import os
import casatools as cc
from os import path

candname = sys.argv[1]
wdir = f"/media/ubuntu/data/localization_processing/{candname}/"

idict = np.load(f"{wdir}/rfc_idict_more.npz",allow_pickle=True)["idict"].item()
for src in idict['jname']:
    try:
        if os.path.exists(f"{wdir}/icrs_image/{src}-image.fits"):
            print(f"Image already made for {src}")
        else:
            myms = glob.glob(f"{wdir}/*{src}.ms")[0]
            applycal(vis=myms,gaintable=[f"{wdir}/{candname}.B0"])
            os.system(f"cd {wdir}/icrs_image; wsclean -weight briggs 0.5 -name {src} -size 2000 2000 -scale 0.3asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam {myms}; cd -")
    except:
        print(f"Could not find {src}")

