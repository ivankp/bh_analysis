#!/bin/bash

source /home/ivanp/.scripts/env.sh

dir=/msu/data/t3work2/ivanp/H2j-ggf
bh=$dir/bh
wt=$dir/wt/8TeV/HT2-unc/CT10nlo
hist=$dir/hist/8TeV/HT2-unc/AntiKt4/CT10nlo/pt30-eta4.4/all

analysis=/home/ivanp/work/bh_analysis
exe=$analysis/bin/hist_Hjets

######

#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$dir/condor/hist/lib
#echo "LD_LIBRARY_PATH=" $LD_LIBRARY_PATH
#which root
#ldd $exe

######

time1=`date +%s`

cp $bh/$1 $TMPDIR/bh_$1
cp $wt/$1 $TMPDIR/wt_$1

time2=`date +%s`

echo "Copying took $((time2-time1)) seconds"

$exe -j2 --bh=$TMPDIR/bh_$1 --wt=$TMPDIR/wt_$1 -o $hist/$1 -c AntiKt4 -s $analysis/config/H3j_13TeV.css

