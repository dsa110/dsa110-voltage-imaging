import numpy as np
from casatasks import *
import sys
import argparse
import glob
import os
import casatools as cc

def bp_gain(bpcal,candname,gaincalib,frb,model_field=False):

    wdir = f"/media/ubuntu/ssd/localization_processing/{candname}/"
    if gaincalib is None:
        gaincalib = glob.glob(wdir+f"/*_{candname}.ms")[0]
    if frb is None:
        frb = wdir+f"/{candname}.ms"

    # run bandpass
    print(f"Bandpass calibration on {bpcal}")
    flagdata(vis=bpcal,mode='unflag')
    bandpass(vis=bpcal,caltable=f"{wdir}/{candname}.B0",refant="104",solint="inf",bandtype="B",combine="obs,scan",uvrange=">0.3klambda",smodel=[7.6,0.,0.,0.])
    
    # apply bandpass 
    print(f"Gain calibration on {gaincalib}")
    applycal(vis=gaincalib,gaintable=[f"{wdir}/{candname}.B0"])
    os.system(f"rm -rf {wdir}/field.ms*")
    os.system(f"rm -rf {wdir}/frb.ms*")
    os.system(f"rm -rf {wdir}/frb_bp.ms*")
    split(vis=gaincalib,outputvis=f"{wdir}/field.ms")

    if model_field:

        # copy model from gaincalib and run gaincal
        a = cc.componentlist()
        a.addcomponent(shape="Point",dir="J2000 01h00m00.00s +70d00m00.00s",flux=1.0,fluxunit='Jy',freq='1.405GHz')
        os.system("rm -rf tmp.complist")
        a.rename('tmp.complist')
        a.close()
        ft(vis=f"{wdir}/field.ms",complist='tmp.complist',spw='0',usescratch=True)
        tb = cc.table()
        tb.open(gaincalib)
        mdl = tb.getcol('MODEL_DATA')
        tb.close()
        tb.open(f"{wdir}/field.ms",nomodify=False)
        tb.putcol('MODEL_DATA',mdl)
        tb.close()
        gaincal(vis=f"{wdir}/field.ms",caltable=f"{wdir}/{candname}.G0",refant="104",gaintype="G",calmode="p",solint="180s",minsnr=3.0,combine="obs,scan",uvrange=">0.3klambda")

    # apply gaincal to field and frb
    print(f"Applying calibration to field.ms and frb.ms")
    if model_field:
        applycal(vis=f"{wdir}/field.ms",gaintable=[f"{wdir}/{candname}.G0"])
        split(vis=frb,outputvis=f"{wdir}/frb.ms",datacolumn="DATA")
        applycal(vis=f"{wdir}/frb.ms",gaintable=[f"{wdir}/{candname}.B0",f"{wdir}/{candname}.G0"],interp=['nearest','linear'])
    split(vis=frb,outputvis=f"{wdir}/frb_bp.ms",datacolumn="DATA")        
    applycal(vis=f"{wdir}/frb_bp.ms",gaintable=[f"{wdir}/{candname}.B0"],interp=['nearest'])
    
def parse_commandline_arguments() -> 'argparse.Namespace':
    """ Parse command line args. """
    parser = argparse.ArgumentParser(
        description="Calibrate FRB.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(
        '--candname',
        type=str,
        nargs='?',        
        help="candname")
    parser.add_argument(
        '--bpcal',
        type=str,
        nargs='?',
        help="bandpass cal ms")
    parser.add_argument(
        '--gaincal',
        type=str,
        nargs='?',
        default=None,
        help="gain cal ms")
    parser.add_argument(
        '--frb',
        type=str,
        nargs='?',
        default=None,
        help="frb ms")
    parser.add_argument(
        '--model_field',
        action='store_true')
    parser.set_defaults(model_field=False)

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    print("Running main program")
    ARGS = parse_commandline_arguments()
    bp_gain(ARGS.bpcal,ARGS.candname,ARGS.gaincal,ARGS.frb,model_field=ARGS.model_field)



