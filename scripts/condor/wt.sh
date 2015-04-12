#!/bin/bash

source /home/ivanp/.scripts/env.sh

dir=/msu/data/t3work2/ivanp/H1j-ggf
bh=$dir/bh
wt=$dir/wt/13TeV/HT2-unc/CT10nlo

BIN=/home/ivanp/work/bh_analysis/bin

######

time1=`date +%s`

cp $bh/$1 $TMPDIR/bh_$file

time2=`date +%s`

echo "Copying took $((time2-time1)) seconds"

######

$BIN/reweigh --bh=$bh/$1 -o $wt/$1 --pdf=CT10nlo -c $BIN/../config/rew_Ht2.xml

