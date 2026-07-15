#!/bin/bash

#candname="220319aaeb"
#ntint="8"
#startoffset="1907"
#stopoffset="1908"
#dm="110.95"
#bpcal="<weights datestring>"

candname=${1}
ntint=${2}
startoffset=${3}
stopoffset=${4}
dm=${5}
bpcal=${6}

odir="/media/ubuntu/ssd/localization_processing/${candname}"

mkdir -p ${odir}

# create burst MS from voltage data

cp -r /dataz/dsa110/candidates/${candname}/Level2/calibration/*${bpcal}* /home/ubuntu/scratch

python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm ${dm} --weights ${bpcal} ${candname}
mv ${odir}/${candname}.ms ${odir}/${candname}_full.ms
python voltages_to_ms.py --ntint ${ntint} --startoffset ${startoffset} --stopoffset ${stopoffset} --dm ${dm} --weights ${bpcal} ${candname}
python massage_ms.py ${candname}_full.ms ${odir}
python massage_ms.py ${candname}.ms ${odir}

# image frb and field ms
mkdir -p ${odir}/images
cd ${odir}/images
wsclean -weight briggs 0.5 -name frb -size 10000 10000 -scale 2asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam ../${candname}.ms
wsclean -weight briggs 0.5 -name frb_full -size 10000 10000 -scale 2asec -minuv-l 200 -niter 20 -auto-threshold 4.0 -theoretic-beam ../${candname}_full.ms
cd /home/ubuntu/proj/dsa110-shell/dsa110-voltage-imaging/scripts








