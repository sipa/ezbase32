#!/bin/bash

CHARSET="0123456789ABCDEFGHIJKLMNOPQRSTUV"

function rnd() {
  echo $(($(dd if=/dev/urandom bs=1 count=4 status=none | hexdump -e '"%u"') % $1))
}

function char() {
  echo "${CHARSET:$(rnd 32):1}"
}

function rchar() {
  echo "${CHARSET:$((1+$(rnd 31))):1}"
}

NUM=$(cat /proc/cpuinfo | fgrep processor | wc -l)

for N in $(seq 1 $NUM); do
  while true; do
    RND="$(char)$(char)$(char)$(char)$(char)$(char)$(char)$(char)$(char)$(char)$(char)$(rchar)"
    echo "$RND: begin at $(date)"
    ./crccollide $RND 8 65
  done &
done

wait
