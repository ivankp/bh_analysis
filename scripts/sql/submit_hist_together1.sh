#!/bin/bash

bin=/home/ivanp/work/bh_analysis/bin
dir=/msu/data/t3work2/ivanp
out=$dir/out

mkdir -p $out

db="$dir/ntuples.db"

sql0="from bh"

for jalg     in "AntiKt4"; do
for dset     in `sqlite3 $db "select distinct dset     from bh"`; do
for particle in `sqlite3 $db "select distinct particle from bh"`; do
for njets    in `sqlite3 $db "select distinct njets    from bh"`; do
for energy   in `sqlite3 $db "select distinct energy   from bh"`; do
for part     in `sqlite3 $db "select distinct part     from bh"`; do
for scales   in `sqlite3 $db "select distinct scales   from wt"`; do


echo $dset $particle $njets $energy $part $scales
echo $db "select dir, file from wt where scales='$scales' and bh_id in (
  select id from bh where dset = '$dset'
    and dset     = '$dset'
    and particle = '$particle'
    and nj       = '$nj'
    and ene      = '$ene'
    and part     = '$part'
"
echo ""


done
done
done
done
done
done
done



echo 'DONE!'

