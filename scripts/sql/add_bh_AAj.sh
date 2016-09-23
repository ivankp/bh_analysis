#!/bin/bash

db=$1

all="/msu/data/t3fast2/giorgi"

set=`sqlite3 $db "select distinct dset from bh" | tail -1`
if [ -z "$set" ]; then
  set=0
fi

for d in `ls $all | grep '^AA[0-9]j$'`; do
  AA=`ls $all/$d`
  dir="$all/$d/$AA"
  echo "$dir"

  arr=(`echo $AA | sed -n \
    's/^AA\([0-9]\)j-\(.*\)-\(.*\)tev-eA\(.*\)-ptA\(.*\)-ej\(.*\)-ptj\(.*\)$/\1 \2 \3 \4 \5 \6 \7/p'`\
  )

  cd $dir
  dir=`pwd -P`
  cd -

  ((set++))

  for p in `ls $dir`; do
    for f in `ls $dir/$p`; do
      echo "$AA/$p/$f"

      sid=`echo $f | sed -n 's/^event\([0-9]*\)\.root$/\1/p'`

      part=`echo $p | tr '_' '-'`

      sqlite3 $db "insert into bh
        (file,dir,dset,particle,njets,energy,part,jet_pt_cut,jet_eta_cut,other_cuts,sid)
        values
        ('$f','$dir/$p',"$set",'AA','${arr[0]}','${arr[2]}','$part','${arr[6]}','${arr[5]}',
         'ptA${arr[4]} etaA${arr[3]}','$sid');"

    done
  done
done
