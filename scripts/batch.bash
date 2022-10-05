#!/bin/bash
#

#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 686.55 220726aabn
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 413.0 --weights 2022-07-30T13:57:33 220801aabd 
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 651.2 220825aaad
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 1146.25 220831aaaj
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 631.05 220914aabz
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 315.0 220920aacl
#python voltages_to_ms.py --ntint 4096 --startoffset 0 --stopoffset 15 --dm 441.53 220926aaeu

python create_ms.py 220825aaad
python create_ms.py 220831aaaj
python create_ms.py 220914aabz
python create_ms.py 220920aacl
python create_ms.py 220926aaeu
