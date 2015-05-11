#!/bin/bash

db=$1

for wt in `sqlite3 -separator '/' $db "select dir, file from wt"`
do
  if [ ! -f "$wt" ]; then
#    date=`ls -l $wt | tr -s ' ' | cut -d' ' -f6,7`
#    echo "$date $wt"
#  else
    echo "No file $wt"
  fi

done

