#!/bin/bash

db=$1

for nj  in 0 1 2 3; do
for ene in 8; do

  for bh in `sqlite3 $db "select id from bh where \
    particle = 'AA' and njets = $nj and energy = $ene"`
  do
    for wt in `sqlite3 $db "select id from wt where bh_id = $bh"`
    do
      arr=(`sqlite3 -separator ' ' $db "select dir, file, scales, pdf from wt where id = $wt"`)

      old="/msu/data/t3work2/ivanp/AA"$nj"j/wt/"$ene"TeV/"${arr[2]}"/"${arr[3]}"/"
      old="$old`sqlite3 $db "select part from bh where id = $bh" | tr '-' '_'`_"
      old="$old`sqlite3 $db "select file from bh where id = $bh"`"

      new="${arr[0]}/${arr[1]}"

      mv -v $old $new
    done
  done

done
done
