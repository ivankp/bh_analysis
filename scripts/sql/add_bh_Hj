#!/bin/bash

db=$1

all="/msu/data/t3work4/luisonig"

scales="HT2"

nsets=`sqlite3 $db "select distinct dset from bh" | tail -1`
if [ -z "$nsets" ]; then
  nsets=0
fi
set=$nsets
proc=""

for d in `ls $all | grep '^H[0-9]jets_ggf$'`; do

  dir="$all/$d/NTuplesFiles"

  cd $dir
  dir=`pwd -P`
  cd -

  proc=""

  for f in `ls $dir`; do
    echo $f
    arr=(`echo $f | sed -n \
's/^H\([0-9]\).0j_\(.*\)_\(.*\)_\([0-9]*\)_pt\([0-9\.]*\)_eta\([0-9\.]*\)_r.*_\(.*\)\.root$/\1 \2 \3 \4 \5 \6 \7/p'`\
    )

    if [ "${arr[1]}" == "GGFHTLS" ]; then continue; fi # skip ntuples with low pt cut

    energy=$((${arr[3]} * 2 / 1000))

    if [ "$proc" != "${arr[1]}" ]; then
      set=`sqlite3 $db "select distinct dset from bh
           where particle = 'H' and njets = ${arr[0]} and process = '${arr[1]}'"`
      echo "set = $set"
      if [ -z "$set" ]; then
        ((nsets++))
        set="$nsets"
      fi
      echo "set = $set"
    fi
    proc="${arr[1]}"

    sqlite3 $db "insert into bh
      (file,dir,dset,particle,njets,process,energy,part,scales,jet_pt_cut,jet_eta_cut,sid)
      values
      ('$f','$dir','$set','H','${arr[0]}','$proc','$energy','${arr[2]}','$scales','${arr[4]}',
       '${arr[5]}','${arr[6]}');"

  done
done

