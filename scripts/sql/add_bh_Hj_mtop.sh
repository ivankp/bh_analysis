#!/bin/bash

db="/msu/data/t3work2/ivanp/ntuples.db"

set=`sqlite3 $db "select distinct dset from bh" | tail -1`
if [ -z "$set" ]; then
  set=0
fi

for j in 1 2 3; do
((set += 1))

dir="/msu/data/t3work4/luisonig/H${j}jets_ggf_mtop/NTuplesFiles"

for f in `find $dir -type f`; do

  file=`basename $f`

  fields=`sed -n "s|RWGT\([0-9]*\)_\([^_]*\)\([0-9]\)\.0j_\([^_]*\)_\([^_]*\)_\([0-9\.]*\)_pt\([0-9\.]*\)_eta\([0-9\.]*\)_r100_\([0-9]*\)\.root|'\2',\3,'\4','\5',\7,\8,\6,\1,\9|p" <<< $file`

  energy=`echo $fields | cut -d',' -f7 | sed 's|.*|scale=1; &/500|' | bc -l`
  sid=`echo $fields | cut -d',' -f8,9 | sed -r 's|(.*),(.*)|scale=0; \1*1000+\2|' | bc -l`

  fields=`echo $fields | cut -d',' -f1-6`
  fields="$fields,$energy,$sid"

  echo $fields

  sqlite3 $db "insert into bh
    (file,dir,dset,particle,njets,process,part,jet_pt_cut,jet_eta_cut,energy,sid)
    values ('$file','$dir',$set,$fields);"

done

done

