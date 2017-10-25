#!/bin/bash

NUM=$(cat /proc/cpuinfo | fgrep processor | wc -l)

for N in $(seq 1 $NUM); do
  while true; do
    ./crccollide "" 8 100
  done &
done | tee -a log_deg13len100.log

wait
