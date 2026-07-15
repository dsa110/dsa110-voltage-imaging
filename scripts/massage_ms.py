from casatasks import *
import sys
import os

ms = sys.argv[1]
odir = sys.argv[2]

split(vis=f"{odir}/{ms}",outputvis=f"{odir}/frb.ms",datacolumn="DATA")
os.system(f"rm -rf {odir}/{ms}; mv {odir}/frb.ms {odir}/{ms}")

