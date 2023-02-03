#!/bin/bash

#candname="220319aaeb"
#ntint="8"
#startoffset="1907"
#stopoffset="1908"
#dm="110.95"
#bpcal="/dataz/dsa110/candidates/220319aaeb/Level2/calibration/2022-03-19_J145907+714019.ms"

candname=${1}
ntint=${2}
startoffset=${3}
stopoffset=${4}
dm=${5}
bpcal=${6}

odir="/media/ubuntu/ssd/localization_processing/${candname}"
mkdir -p ${odir}

# create burst MS from voltage data
python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm ${dm} ${candname}
mv ${odir}/${candname}.ms ${odir}/${candname}_full.ms
python voltages_to_ms.py --ntint ${ntint} --startoffset ${startoffset} --stopoffset ${stopoffset} --dm ${dm} ${candname}
# create field MS from realtime visibilities
python create_ms.py ${candname}
# bandpass calibration of burst and field MSs
python bp_gain_cal.py --candname ${candname} --bpcal ${bpcal}
# create MSs on VLBI calibrators before and after the burst. 
python create_all_icrs.py ${candname}
# calibrate and image VLBI MSs
python calibrate_image_icrs.py ${candname}
# image frb and field ms
mkdir -p ${odir}/images
cd ${odir}/images
wsclean -weight briggs 0.5 -name frb -size 10000 10000 -scale 2asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam ../frb_bp.ms
wsclean -weight briggs 0.5 -name fimg -size 6000 6000 -scale 2asec -minuv-l 200 -niter 1000 -auto-threshold 4.0 -theoretic-beam ../field.ms
cd /home/ubuntu/proj/dsa110-shell/dsa110-voltage-imaging/scripts








