#!/bin/bash

source /home/ivanp/.scripts/env.sh

file=$1

nj=0
if [[ $file == "H1.0j_"* ]]; then
  nj=1
elif [[ $file == "H2.0j_"* ]]; then
  nj=2
elif [[ $file == "H3.0j_"* ]]; then
  nj=3
fi

ene=""
if [[ $file == *"_4000_"* ]]; then
  ene=8TeV
elif [[ $file == *"_6500_"* ]]; then
  ene=13TeV
fi

scales=HT2-unc
pdf=CT10nlo
jalg=AntiKt4
cuts=pt30-eta4.4

dir=/msu/data/t3work2/ivanp/H"$nj"j-ggf
bh=$dir/bh
wt=$dir/wt/$ene/$scales/$pdf
hist=$dir/hist/$ene/$scales/$jalg/$pdf/$cuts/all

analysis=/home/ivanp/work/bh_analysis
exe=$analysis/bin/hist_Hjets

######

#LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$dir/condor/hist/lib
#echo "LD_LIBRARY_PATH=" $LD_LIBRARY_PATH
#which root
#ldd $exe

######

time1=`date +%s`

cp $bh/$1 $TMPDIR/bh_$file
cp $wt/$1 $TMPDIR/wt_$file

time2=`date +%s`

echo "Copying took $((time2-time1)) seconds"

$exe -j $nj --bh=$TMPDIR/bh_$file --wt=$TMPDIR/wt_$file -o $hist/$file -c $jalg

