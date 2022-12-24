#!/bin/bash
#

#python create_all_icrs.py 221116aaab 01h24m50.63s +72d39m13.5s
#python create_all_icrs.py 221113aaao 04h45m38.64s +70d18m26.6s
#python create_all_icrs.py 221101aaar 05h04m14.77s +72d26m49.0s
#python create_all_icrs.py 221101aaag 22h48m51.79s +70d40m53.4s
#python create_all_icrs.py 221029aado 09h27m51.24s +72d27m09.4s
#python create_all_icrs.py 221027aags 08h43m29.28s +72d06m03.4s
#python create_all_icrs.py 221012aaab 18h43m12.1s +70d31m27.6s

python do3.py 221116aaab
python do3.py 221113aaao
python do3.py 221101aaar
python do3.py 221101aaag
python do3.py 221029aado
python do3.py 221027aags

# create burst MS from voltage data, using optimal integration, sample, and DM
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 602.25 221203aaaa
# create field MS from realtime visibilities
#python create_ms.py 221203aaaa
# bandpass and gain calibration of burst and field MSs
# typically the bandpass cal is the most recent 1459+716 transit in
# /dataz/dsa110/operations/calibration/*1459+716*.ms
#python bp_gain_cal.py --candname 221203aaaa --bpcal /dataz/dsa110/operations/calibration/2022-12-03_1459+716.ms
# create MSs on VLBI calibrators before and after the burst. 
#python create_icrs.py 220914aabz <RA> <DEC>




