#!/bin/bash

#candname="220319aaeb"
#ntint="8"
#startoffset="1907"
#stopoffset="1908"
#dm="110.95"
#bpcal="/dataz/dsa110/candidates/220319aaeb/Level2/calibration/2022-03-19_J145907+714019.ms"

candname=${1}

odir="/media/ubuntu/data/localization_processing/${candname}"

# create burst MS from voltage data
python calibrate_image_icrs.py ${candname}
# image frb and field ms
cd ${odir}/images
wsclean -weight briggs 0.5 -name frb -size 10000 10000 -scale 2asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam ../frb_bp.ms
wsclean -weight briggs 0.5 -name fimg -size 6000 6000 -scale 2asec -minuv-l 200 -niter 1000 -auto-threshold 4.0 -theoretic-beam ../field.ms
cd /home/ubuntu/proj/dsa110-shell/dsa110-voltage-imaging/scripts








