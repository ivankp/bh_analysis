#!/bin/bash

bin=/home/ivanp/work/bh_analysis/bin
dir=/msu/data/t3work2/ivanp
out=$dir/out

mkdir -p $out

jalg="AntiKt4"

db="$dir/ntuples.db"

sql="select bh_id, id from wt"
if [[ $1 ]]; then
  sql="$sql $1"
fi

for ids in `sqlite3 $db "$sql"`
do

  bh=`echo $ids | cut -d'|' -f1`
  wt=`echo $ids | cut -d'|' -f2`

  bh=(`sqlite3 -separator ' ' "$db" "select dir, file, dset, particle, njets from bh where id = $bh"`)
  wt=(`sqlite3 -separator ' ' "$db" "select dir, file from wt where id = $wt"`)

  base="`basename ${wt[1]} .root`_$jalg"

  hist="$dir/hist/${bh[2]}_${bh[3]}${bh[4]}j/$base.root"

  if [[ -f "$hist" ]]; then continue; fi

  echo $hist

  if [[ "${bh[3]}" == "H" ]]; then
    exe=hist_Hjets
  elif [[ "${bh[3]}" == "AA" ]]; then
    exe=hist_AAjets
  else
    echo "Unexpected particle ${bh[3]}"
    continue
  fi

  echo "
Universe   = vanilla
Executable = $bin/$exe
Error      = $out/$base.err
Output     = $out/$base.out
Log        = $out/$base.log

getenv = True
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = ${bh[0]}/${bh[1]}, ${wt[0]}/${wt[1]}

Arguments = -j ${bh[4]} --bh=${bh[0]}/${bh[1]} --wt=${wt[0]}/${wt[1]} -o $hist -c $jalg
Queue
" | condor_submit - > /dev/null

#sleep 5

done

echo 'DONE!'

