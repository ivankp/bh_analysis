#!/usr/bin/env python

import sys, sqlite3

if len(sys.argv) != 2 :
  print 'usage: %s ntuples.db' % sys.argv[0]
  sys.exit(1)

_dir='/msu/data/t3work2/ivanp/wt'

db = sqlite3.connect(sys.argv[1])
cur = db.cursor()

cur.execute("select id, particle, njets, energy, part, sid, dset FROM bh")

pnj = ""

for scales in ['HT2-unc']:
  for pdf in ['CT10nlo']:

    for bh in cur.fetchall():
      f = '%s%dj_%dTeV_%s_%s_%s_%s.root' % ( bh[1:-1] + (scales,pdf) )

      _pnj = '%s%dj' % bh[1:3]
      if pnj != _pnj:
        pnj = _pnj
        print pnj

      cur.execute(
        "insert into wt (bh_id,file,dir,scales,pdf) values ('%s','%s','%s/%d_%s%dj','%s','%s')" % (
          bh[0], f, _dir, bh[-1], bh[1], bh[2], scales, pdf )
      )

db.commit()

