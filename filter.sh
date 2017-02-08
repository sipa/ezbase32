#!/bin/bash

CONST="1-5,1-15,0 1-6,1-7,0 1-4,1-71,0 1-5,1-71,0.8912 1-5,1-59,0.8616 1-5,1-39,0.8118 1-8,1-71,1.073"

/usr/bin/time ./crcmiddle 5 15 15 $CONST | cut -d ' ' -f 1 | uniq | tee phase_1_"$1".log | buffer |
/usr/bin/time ./crcmiddle 6  7  7 $CONST | cut -d ' ' -f 1 | uniq | tee phase_2_"$1".log | buffer |
/usr/bin/time ./crcmiddle 4 50 50 $CONST | cut -d ' ' -f 1 | uniq | tee phase_3_"$1".log | buffer |
/usr/bin/time ./crcmiddle 6 15 15 $CONST | cut -d ' ' -f 1 | uniq | tee phase_4_"$1".log | buffer |
/usr/bin/time ./crcmiddle 5 39 39 $CONST | cut -d ' ' -f 1 | uniq | tee phase_5_"$1".log | buffer |
/usr/bin/time ./crcmiddle 4 71 71 $CONST | cut -d ' ' -f 1 | uniq | tee phase_6_"$1".log | buffer |
/usr/bin/time ./crcmiddle 5 71 71 $CONST | cut -d ' ' -f 1 | uniq | tee phase_7_"$1".log | buffer |
/usr/bin/time ./crcmiddle 6 71 71 $CONST | cut -d ' ' -f 1 | uniq | tee phase_8_"$1".log | buffer
/usr/bin/time ./crcmiddle 8 10 10 $CONST | cut -d ' ' -f 1 | uniq >phase_9_"$1".log
