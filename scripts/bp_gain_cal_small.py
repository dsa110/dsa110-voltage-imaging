import numpy as np
from casatasks import *
import sys
import argparse
import glob
import os
import casatools as cc

def bp(bpcal,candname,flux_Jy=7.6):

    wdir = f"/media/ubuntu/ssd/localization_processing/{candname}/"

    # run bandpass
    print(f"Bandpass calibration on {bpcal}")
    flagdata(vis=bpcal,mode='unflag')
    bandpass(vis=bpcal,caltable=f"{wdir}/{candname}.B0",refant="103",solint="inf",bandtype="B",combine="obs,scan",uvrange=">0.3klambda",smodel=[flux_Jy,0.,0.,0.])
    
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

    args = parser.parse_args()
    return args

if __name__ == '__main__':
    print("Running main program")
    ARGS = parse_commandline_arguments()
    bp(ARGS.bpcal,ARGS.candname)



