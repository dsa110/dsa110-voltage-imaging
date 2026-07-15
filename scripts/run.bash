#!/bin/bash
#

#source batch.bash 220204aaai 32 476 477 612.2 /dataz/dsa110/candidates/220204aaai/Level2/calibration/2022-02-04_J145907+714019.ms
#source batch.bash 220207aabh 16 953 959 262.3 /dataz/dsa110/candidates/220207aabh/Level2/calibration/2022-02-07_J145907+714019.ms
#source batch.bash 220208aaaa 64 238 239 437.0 /dataz/dsa110/candidates/220208aaaa/Level2/calibration/2022-02-09_J145907+714019.ms
#source batch.bash 220307aaae 8 1906 1907 499.15 /dataz/dsa110/candidates/220307aaae/Level2/calibration/2022-03-08_J145907+714019.ms
#source batch.bash 220310aaam 32 476 478 462.15 /dataz/dsa110/candidates/220310aaam/Level2/calibration/2022-03-09_J145907+714019.ms
#source batch.bash 220330aaan 256 59 60 468.1 /dataz/dsa110/candidates/220330aaan/Level2/calibration/2022-03-30_J145907+714019.ms
#source batch.bash 220418aaai 8 1907 1908 623.45 /dataz/dsa110/candidates/220418aaai/Level2/calibration/2022-04-18_J145907+714019.ms
#source batch.bash 220424aabq 4 3815 3817 863.35 /dataz/dsa110/candidates/220424aabq/Level2/calibration/2022-04-24_J145907+714019.ms
#source batch.bash 220506aabd 16 952 953 396.93 /dataz/dsa110/candidates/220506aabd/Level2/calibration/2022-05-06_J145907+714019.ms
#source batch.bash 220825aaad 8 1905 1906 651.2 /dataz/dsa110/operations/calibration/2022-08-25_1459+716.ms
#source batch.bash 220831aaaj 64 238 239 1146.25 /dataz/dsa110/operations/calibration/2022-08-31_1459+716.ms
#source batch.bash 220914aabz 4 3814 3815 631.05 /dataz/dsa110/operations/calibration/2022-09-14_1459+716.ms
#source batch.bash 220920aacl 16 953 954 315.0 /dataz/dsa110/operations/calibration/2022-09-20_1459+716.ms
#source batch.bash 220926aaeu 4 3815 3817 441.53 /dataz/dsa110/operations/calibration/2022-09-25_1459+716.ms
#source batch.bash 221002aaab 512 29 30 322.5 /dataz/dsa110/operations/calibration/2022-10-01_1459+716.ms
#source batch.bash 221012aaab 32 476 477 442.2 /dataz/dsa110/operations/calibration/2022-10-12_1459+716.ms
#source batch.bash 221027aags 8 1906 1907 452.5 /dataz/dsa110/operations/calibration/2022-10-27_1459+716.ms
#source batch.bash 221029aado 32 477 478 1391.05 /dataz/dsa110/operations/calibration/2022-10-28_1459+716.ms
#source batch.bash 221101aaag 4 3807 3810 490.7 /dataz/dsa110/operations/calibration/2022-11-01_1459+716.ms
#source batch.bash 221101aaar 16 951 954 1476.5 /dataz/dsa110/operations/calibration/2022-11-01_1459+716.ms
#source batch.bash 221113aaao 4 3813 3814 411.4 /dataz/dsa110/operations/calibration/2022-11-10_1459+716.ms
#source batch.bash 221116aaab 64 238 239 640.6 /dataz/dsa110/operations/calibration/2022-11-16_1459+716.ms
#source batch.bash 221203aaaa 8 1906 1908 602.25 /dataz/dsa110/operations/calibration/2022-12-03_1459+716.ms
#source batch.bash 221219aabz 256 59 63 706.7 /dataz/dsa110/operations/calibration/2022-12-18_1459+716.ms
#source batch.bash 221224aaab 64 237 240 745.3 /dataz/dsa110/operations/calibration/2022-12-24_1459+716.ms
#source batch.bash 220426aaaw 8 1907 1908 269.5 /dataz/dsa110/candidates/220426aaaw/Level2/calibration/2022-05-08_J145907+714019.ms

#source batch2.bash 220726aabn 32 476 477 686.55 /dataz/dsa110/candidates/220726aabn/Level2/calibration/220726aabn.B0
#source batch.bash 221018aaaj 256 59 60 219.46 /dataz/dsa110/operations/calibration/2022-10-17_1459+716.ms
#source batch.bash 221025aanu 256 59 60 219.46 /dataz/dsa110/operations/calibration/2022-10-24_1459+716.ms

#source batch.bash 230112aaaw 64 237 239 650.2 /dataz/dsa110/operations/calibration/2023-01-12_1459+716.ms
#source batch.bash 230124aaav 128 118 120 590.0 /dataz/dsa110/operations/calibration/2023-01-24_1459+716.ms
#source batch.bash 230202aada 32 475 477 810.0 /dataz/dsa110/operations/calibration/2023-02-02_1229+020.ms
#source batch.bash 230216aaag 32 476 477 828.0 /dataz/dsa110/operations/calibration/2023-02-16_1229+020.ms
#source batch.bash 230307aaao 16 952 956 610.15 /dataz/dsa110/operations/calibration/2023-03-06_1459+716.ms
#source batch.bash 230325aaag 4 3818 3821 912.38 /dataz/dsa110/operations/calibration/2023-03-25_1459+716.ms
#source batch.bash 230501aaae 8 1905 1907 533.7 /dataz/dsa110/operations/calibration/2023-05-01_1459+716.ms
#source batch.bash 230521aaad 16 952 953 1345.65 /dataz/dsa110/operations/calibration/2023-05-22_1459+716.ms
#source batch.bash 230626aaai 8 1908 1909 451.2 /dataz/dsa110/operations/calibration/2023-06-26_1459+716.ms
#source batch.bash 230628aaaa 32 476 477 345.15 /dataz/dsa110/operations/calibration/2023-06-26_1459+716.ms
#source batch.bash 230712aadj 2 7627 7630 586.96 /dataz/dsa110/operations/calibration/2023-07-12_1459+716.ms
#source batch.bash 230714aacu 4 3815 3817 838.95 /dataz/dsa110/operations/calibration/2023-07-14_1459+716.ms
#source batch.bash 230814aaas 16 952 954 696.35 /dataz/dsa110/operations/calibration/2023-08-14_1459+716.ms
#source batch.bash 230913aaao 16 953 954 518.7 /dataz/dsa110/operations/calibration/2023-09-12_1459+716.ms
#source batch.bash 230928aaag 16 237 239 264.5 /dataz/dsa110/operations/calibration/2023-09-27_1459+716.ms
#source batch.bash 231120aabi 128 118 121 438.9 /dataz/dsa110/operations/calibration/2023-11-20_1459+716.ms
#source batch.bash 231123aacm 32 476 477 396.7 /dataz/dsa110/operations/calibration/2023-11-24_1459+716.ms
#source batch.bash 231129aaaj 32 476 477 720.6 /dataz/dsa110/operations/calibration/2023-11-28_1459+716.ms
#source batch.bash 231210aaaa 4 3809 3812 1131.65 /dataz/dsa110/operations/calibration/2023-12-09_1459+716.ms
#source batch.bash 231215aaat 8 1903 1905 977.05 /dataz/dsa110/operations/calibration/2023-12-13_1459+716.ms
#source batch.bash 231220aaae 16 953 954 491.2 /dataz/dsa110/operations/calibration/2023-12-19_1459+716.ms
#source batch.bash 240104aaab 64 237 239 1351.0 /dataz/dsa110/operations/calibration/2024-01-03_1459+716.ms
#source batch.bash 240123aaaf 256 59 60 1462.0 /dataz/dsa110/operations/calibration/2024-01-23_1459+716.ms
#source batch.bash 240119aacg 64 238 240 483.1  /dataz/dsa110/operations/calibration/2024-01-19_1459+716.ms
#source batch.bash 240122aaag 4 3806 3808 959.9  /dataz/dsa110/operations/calibration/2024-01-22_1459+716.ms
#source batch.bash 240203aacl 16 953 954 272.25 /dataz/dsa110/operations/calibration/2024-02-03_1459+716.ms
#source batch.bash 240213aani 32 476 479 357.4 /dataz/dsa110/operations/calibration/2024-02-13_1459+716.ms
#source batch.bash 240215aawc 32 476 478 549.5 /dataz/dsa110/operations/calibration/2024-02-15_1459+716.ms
#source batch.bash 240224aaad 4 3810 3812 881.2 /dataz/dsa110/operations/calibration/2024-02-23_1459+716.ms
#source batch.bash 240229aaad 8 1907 1908 491.15 /dataz/dsa110/operations/calibration/2024-02-29_1459+716.ms
source batch.bash 240304aaax 64 236 238 443.2 /dataz/dsa110/operations/calibration/2024-03-04_1459+716.ms










