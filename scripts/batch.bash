#!/bin/bash

candname="220319aaeb"
ntint="8"
startoffset="1907"
stopoffset="1908"
dm="110.95"
bpcal="/dataz/dsa110/candidates/220319aaeb/Level2/calibration/2022-03-19_J145907+714019.ms"

odir="/media/ubuntu/ssd/localization_processing/${candname}"

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






