import numpy as np
from casatasks import *
import sys
import glob
import os
import casatools as cc

candname = sys.argv[1]
wdir = f"/media/ubuntu/ssd/localization_processing/{candname}/"
os.system(f"mkdir -p {wdir}/icrs_image")

idict = np.load(f"{wdir}/rfc_idict.npz",allow_pickle=True)["idict"].item()
for src in idict['jname']:

    myms = glob.glob(f"{wdir}/*{src}.ms")[0]
    applycal(vis=myms,gaintable=[f"{wdir}/{candname}.B0"])
    os.system(f"cd {wdir}/icrs_image; wsclean -weight briggs 0.5 -name {src} -size 2000 2000 -scale 0.3asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam {myms}; cd -")

