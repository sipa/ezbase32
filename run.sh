#!/bin/bash

NUM=$(cat /proc/cpuinfo | fgrep processor | wc -l)

for N in $(seq 1 $NUM); do
  while true; do
    ./crccollide "" 8 105
  done &
done

wait
