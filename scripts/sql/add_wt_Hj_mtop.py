#!/usr/bin/env python

import sqlite3

_dir='/msu/data/t3work2/ivanp/wt'

db = sqlite3.connect('/msu/data/t3work2/ivanp/ntuples.db')
cur = db.cursor()

cur.execute("select id, particle, njets, energy, part, sid, dset FROM bh where dset = 10 or dset = 11 or dset = 12")

pnj = ""

for scales in ['HT2-unc']:
  for pdf in ['CT10nlo']:

    for bh in cur.fetchall():
      f = '%s%dj_%dTeV_%s_%s_%s_%s.root' % ( bh[1:-1] + (scales,pdf) )

      _pnj = '%s%dj' % bh[1:3]
      if pnj != _pnj:
        pnj = _pnj
        print pnj

      sql = "insert into wt (bh_id,file,dir,scales,pdf) values ('%s','%s','%s/%d_%s%dj_mtop','%s','%s')" % (
          bh[0], f, _dir, bh[-1], bh[1], bh[2], scales, pdf )

      cur.execute(sql)

db.commit()

