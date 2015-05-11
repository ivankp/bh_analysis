#!/bin/bash

db=$1

if [ ! -f $db ]; then
  sqlite3 $db "create table bh (
    id INTEGER PRIMARY KEY,
    file TEXT,
    dir TEXT,
    dset INTEGER,
    particle TEXT,
    njets INTEGER,
    process TEXT,
    energy REAL,
    part TEXT,
    scales TEXT,
    pdf TEXT,
    jet_pt_cut REAL,
    jet_eta_cut REAL,
    other_cuts TEXT,
    sid INTEGER
  );"
  sqlite3 $db "create table wt (
    id INTEGER PRIMARY KEY,
    bh_id INTEGER,
    file TEXT,
    dir TEXT,
    scales TEXT,
    pdf TEXT
  );"
#  sqlite3 $db "create table hist (
#    wt_id INTEGER,
#    file TEXT,
#    dir TEXT,
#    jet TEXT
#  );"
  sqlite3 $db ".dump"
fi
