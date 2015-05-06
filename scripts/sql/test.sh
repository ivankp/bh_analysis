#!/bin/bash

sqlite3 ntuples.db "insert into bh
        (file,dir,particle,njets,energy,part,jet_pt_cut,jet_eta_cut,other_cuts)
        values
        ('event197.root','/msu/data/t3fast2/giorgi/AA2j/AA2j-v2.1.1-8tev-eA3-ptA20-ej5-ptj20/RS_1','AA','2','8','RS_1','20','5','ptA20 etaA3')
      ;"
