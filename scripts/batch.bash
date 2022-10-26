#!/bin/bash
#

python voltages_to_ms.py --ntint 256 --startoffset 59 --stopoffset 60 --dm 235.0 221025aanu
python create_ms.py 221025aanu
python bp_gain_cal.py --candname 221025aanu --bpcal /dataz/dsa110/operations/calibration/2022-10-24_1459+716.ms
