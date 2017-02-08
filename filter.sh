#!/bin/bash

CONST="1-5,1-15,0 1-6,1-7,0 1-4,1-71,0 1-5,1-71,0.8912 1-5,1-59,0.8616 1-5,1-39,0.8118 1-8,1-71,1.073"

./crcmiddle $2 5 15 15 $CONST | tee phase_1.log |
./crcmiddle $2 56 6  7  7 $CONST | tee phase_2.log |
./crcmiddle $2 56 4 50 50 $CONST | tee phase_3.log |
./crcmiddle $2 6 15 15 $CONST | tee phase_4.log |
./crcmiddle $2 5 39 39 $CONST | tee phase_5.log |
./crcmiddle $2 4 71 71 $CONST | tee phase_6.log |
./crcmiddle $2 5 71 71 $CONST | tee phase_7.log |
./crcmiddle $2 6 71 71 $CONST | tee phase_8.log |
./crcmiddle $2 8 10 10 $CONST | cut -d ' ' -f 1 | uniq >phase_9.log
