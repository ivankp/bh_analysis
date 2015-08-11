#!/bin/bash

db=~/disk2/ntuples.db

for f in `sqlite3 -separator ',' $db "
  SELECT wt.*
  FROM wt
  JOIN bh ON wt.bh_id = bh.id
  WHERE wt.pdf = 'CT10nlo' and wt.scales = 'HT2-unc' and bh.particle = 'H'
"`; do
  sql="insert into wt (bh_id,file,dir,scales,pdf) values (`echo $f | sed 's/^[^,]*,//;s/,\([^,]*\)/,\"\1\"/g;s/CT10nlo/MMHT2014nlo/g'`)"
  echo "$sql"
  sqlite3 $db "$sql"
done

