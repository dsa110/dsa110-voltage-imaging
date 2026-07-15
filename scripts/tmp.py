import numpy as np
from casatasks import *
import sys
import glob
import os
import casatools as cc

candname = sys.argv[1]
wdir = f"/media/ubuntu/ssd/localization_processing/{candname}/"
wdir_tmp = f"/media/ubuntu/ssd/localization_processing/{candname}_tmp/"

# recalibrate frb data
split(vis=f"{wdir}/frb_bp.ms",outputvis=f"{wdir}/frb_bp_cal.ms")
split(vis=f"{wdir}/field.ms",outputvis=f"{wdir}/field_cal.ms",width=4,datacolumn="data")
applycal(vis=f"{wdir}/frb_bp_cal.ms",gaintable=[f"{wdir_tmp}/2326.G0"])
applycal(vis=f"{wdir}/field_cal.ms",gaintable=[f"{wdir_tmp}/2326.G0"])

# recalibrate and image icrs data
#idict = np.load(f"{wdir}/rfc_idict.npz",allow_pickle=True)["idict"].item()
#for src in idict['jname']:
#    try:
#        myms = glob.glob(f"{wdir}/*{src}.ms")[0]
#        os.system(f"rm -rf {wdir}/tmp.ms*")
#        split(vis=myms,outputvis=f"{wdir}/tmp.ms",width=16)
#        applycal(vis=f"{wdir}/tmp.ms",gaintable=[f"{wdir}/2326.G0"])
#        os.system(f"cd {wdir}/icrs_image; wsclean -weight briggs 0.5 -name {src} -size 2000 2000 -scale 0.3asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam ../tmp.ms; cd -")
        
#    except:
#        print(f"Could not find {src}")

