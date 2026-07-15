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
obsdec=${7}

# source tmp.bash 250924mptq 2 1790 1792 405.25 na 54.5

odir="/media/ubuntu/ssd/localization_processing/${candname}"
mkdir -p ${odir}

# create burst MS from voltage data
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm ${dm} --obsdec ${obsdec} ${candname}
#mv ${odir}/${candname}.ms ${odir}/${candname}_full.ms
#python voltages_to_ms.py --ntint ${ntint} --startoffset ${startoffset} --stopoffset ${stopoffset} --dm ${dm} --obsdec ${obsdec} ${candname}
#mv ${odir}/${candname}.ms ${odir}/${candname}_frb.ms
### create field MS from realtime visibilities
#python create_ms.py ${candname} model
### bandpass calibration of burst and field MSs
#python bp_gain_cal_small.py --candname ${candname} --bpcal ${bpcal}
#python gain_cal.py --candname ${candname} --bpcal ${bpcal}

### create MSs on VLBI calibrators before and after the burst. 
#python create_all_icrs.py ${candname}
### calibrate and image VLBI MSs
#python calibrate_image_icrs.py ${candname}
## image frb and field ms
#mkdir -p ${odir}/images

#cd ${odir}
#wsclean -weight briggs 0.5 -name frb -size 10000 10000 -scale 2asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam frb.ms
#wsclean -weight briggs 0.5 -name fimg -size 6000 6000 -scale 2asec -minuv-l 200 -niter 1000 -auto-threshold 4.0 -theoretic-beam field.ms
#wsclean -weight briggs 0.5 -name frbfield -size 6000 6000 -scale 2asec -minuv-l 200 -niter 1000 -auto-threshold 4.0 -theoretic-beam full.ms
#cd /home/ubuntu/proj/dsa110-shell/dsa110-voltage-imaging/scripts








