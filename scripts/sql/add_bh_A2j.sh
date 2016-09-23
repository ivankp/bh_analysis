#!/bin/bash

db="/msu/data/t3work2/ivanp/ntuples.db"

all="/msu/data/t3work1/Y2j"

set=`sqlite3 $db "select distinct dset from bh" | tail -1`
if [ -z "$set" ]; then
  set=0
fi
((set += 1))

for f in `find $all -type f`; do

  fields=`sed -n "s|^\(.*\)/\(Y\([0-9]\+\)j_\([0-9]\+\)TeV_\([^0-9]\+\)\([0-9]\+\)*_Et\([0-9\.]\+\)GeV_\(.\+\)\.root\)|'\2','\1',$set,'A','\3','\4','\5','\7','\8'|p" <<< $f \
         | sed "s/'R'/'RS'/"`

  sqlite3 $db "insert into bh
    (file,dir,dset,particle,njets,energy,part,jet_eta_cut,sid)
    values ($fields);"

done

