#!/bin/bash

all="/msu/data/t3fast2/giorgi"

for d in `ls $all | grep 'AA[0-9]j'`; do
  AA=`ls $all/$d`
  dir="$all/$d/$AA"
  echo "$dir"
  echo $AA | sed -n 's/^AA\([0-9]\)j-\(.*\)-\(.*\)tev-eA\(.*\)-ptA\(.*\)-ej\(.*\)-ptj\(.*\)$/\1 \2 \3 \4 \5 \6 \7/p'
done

