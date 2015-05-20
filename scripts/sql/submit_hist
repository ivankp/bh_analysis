#!/bin/bash

export dir=/msu/data/t3work2/ivanp
export out=$dir/out

mkdir -p $out
mkdir -p $dir/hist

export db="$dir/ntuples.db"

#####################################################################

for jalg in AntiKt4
do

for VBF in none
do
export VBF

for set in `sqlite3 $db "
  SELECT distinct particle, njets, energy, part, wt.scales
  FROM bh
  JOIN wt ON bh.id = wt.bh_id
  WHERE wt.scales = 'HT2-unc' and particle = 'H'
"`; do

arr=(`echo $set | tr '|' ' '`)
arr[2]=`printf "%.0f" ${arr[2]}` # round real-valued energy
base="${arr[0]}${arr[1]}j_${arr[2]}TeV_${arr[3]}_${arr[4]}_$jalg"

if [ "$VBF" != "none" ]; then
  base=$base"_VBF$VBF"
fi

echo $base

echo "
Universe   = vanilla
Executable = `pwd -P`/hist_Hjets
Error      = $out/$base.err
Output     = $out/$base.out
Log        = $out/$base.log

getenv = True

Arguments = $base
Queue
" | condor_submit - > /dev/null

done
done
done

echo 'DONE!'

