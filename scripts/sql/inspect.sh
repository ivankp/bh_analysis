#!/bin/bash

printf "num rows "
sqlite3 ntuples.db "SELECT Count(*) FROM bh"

for col in particle dset process energy part jet_pt_cut jet_eta_cut other_cuts
do
  echo "*** $col ***"
  sqlite3 ntuples.db "select distinct $col from bh"
done

