#!/bin/bash

sql_find () {
  sqlite3 -separator '/' ntuples.db "
    SELECT dir, file
    FROM bh
    WHERE particle == 'AA'
      AND njets == $1
      AND part == '$2'
      AND sid == $3
  ;"
}

sql_rm () {
  sql_find $1 $2 $3
  sqlite3 -separator '/' ntuples.db "
    DELETE
    FROM bh
    WHERE particle == 'AA'
      AND njets == $1
      AND part == '$2'
      AND sid == $3
  ;"
}

sql_rm 2 RS-1 76
sql_rm 3 RS-1 529
sql_rm 3 RS-1 565
sql_rm 3 RS-2 513
sql_rm 3 RS-2 530
sql_rm 3 RS-2 532
sql_rm 3 RS-2 540
sql_rm 3 RS-2 543
sql_rm 3 RS-2 544
sql_rm 3 RS-2 545
sql_rm 3 RS-2 546
sql_rm 3 RS-2 549
sql_rm 3 RS-2 552
sql_rm 3 RS-2 569
sql_rm 3 RS-2 579
sql_rm 3 RS-2 584
sql_rm 3 RS-2 585
sql_rm 3 RS-2 598
sql_rm 3 RS-2 599
sql_rm 3 RS-2 601
sql_rm 3 RS-2 602
sql_rm 3 RS-2 608
sql_rm 3 RS-2 616
sql_rm 3 RS-2 622
sql_rm 3 RS-2 624
sql_rm 3 RS-2 629
sql_rm 3 RS-2 632
sql_rm 3 RS-2 637
sql_rm 3 RS-2 651
sql_rm 3 RS-2 655
sql_rm 3 RS-2 661
sql_rm 3 RS-2 662
sql_rm 3 RS-2 667
sql_rm 3 RS-2 671
sql_rm 3 RS-2 681
sql_rm 3 RS-2 685
sql_rm 3 RS-2 694
sql_rm 3 RS-2 700
sql_rm 3 RS-2 711
sql_rm 3 RS-2 714
sql_rm 3 RS-2 728
sql_rm 3 RS-2 729
sql_rm 3 RS-2 734
sql_rm 3 RS-2 763
sql_rm 3 RS-2 766
sql_rm 3 RS-2 780
sql_rm 3 RS-2 784
sql_rm 3 RS-2 785
sql_rm 3 RS-2 789
sql_rm 3 RS-2 790
sql_rm 3 RS-2 797
sql_rm 3 RS-2 805
sql_rm 3 RS-2 829
sql_rm 3 RS-2 830
sql_rm 3 RS-2 839
sql_rm 3 RS-2 841
sql_rm 3 RS-2 844
sql_rm 3 RS-2 846
sql_rm 3 RS-2 856
sql_rm 3 RS-2 860
sql_rm 3 RS-2 863
sql_rm 3 RS-2 870
sql_rm 3 RS-2 872
sql_rm 3 RS-2 893
sql_rm 3 RS-2 896
sql_rm 3 RS-2 917
sql_rm 3 RS-2 931
sql_rm 3 RS-2 932
sql_rm 3 RS-2 934
sql_rm 3 RS-2 935
sql_rm 3 RS-2 938
sql_rm 3 RS-2 939
sql_rm 3 RS-2 941
sql_rm 3 RS-2 947
sql_rm 3 RS-2 960
sql_rm 3 RS-2 974
sql_rm 3 RS-2 977
sql_rm 3 RS-2 978
sql_rm 3 RS-2 982
sql_rm 3 RS-2 989
sql_rm 3 RS-2 992
