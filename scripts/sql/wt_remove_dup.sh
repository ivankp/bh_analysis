#!/bin/bash

db=$1

tmp="/tmp/wt_dump.txt"
clean="/tmp/wt_dump_clean.txt"

if [ ! -f $db ]; then

  sqlite3 $db "select bh_id,file,dir,scales,pdf from wt" > $tmp

  sqlite3 $db "drop table wt"
  sqlite3 $db "vacuum"

  sqlite3 $db "create table wt (
    id INTEGER PRIMARY KEY,
    bh_id INTEGER,
    file TEXT,
    dir TEXT,
    scales TEXT,
    pdf TEXT
  );"

fi

awk '!seen[$0]++' $tmp | nl -w1 -s'|' > $clean

# sqlite> .separator ","
# sqlite> .import /tmp/wt_dump_clean.txt test
