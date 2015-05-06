#!/bin/bash

sql_find () {
  sqlite3 -separator '/' ntuples.db "
    SELECT dir, file
    FROM bh
    WHERE particle IS 'AA'
      AND njets == $1
      AND part IS '$2'
      AND file IS 'event$3.root'
  ;"
}

sql_rm () {
  sqlite3 -separator '/' ntuples.db "
    DELETE
    FROM bh
    WHERE particle IS 'AA'
      AND njets == $1
      AND part IS '$2'
      AND file IS 'event$3.root'
  ;"
}

sql_rm 2 RS_1 76
sql_rm 3 RS_1 529
sql_rm 3 RS_1 565
sql_rm 3 RS_2 513
sql_rm 3 RS_2 530
sql_rm 3 RS_2 532
sql_rm 3 RS_2 540
sql_rm 3 RS_2 543
sql_rm 3 RS_2 544
sql_rm 3 RS_2 545
sql_rm 3 RS_2 546
sql_rm 3 RS_2 549
sql_rm 3 RS_2 552
sql_rm 3 RS_2 569
sql_rm 3 RS_2 579
sql_rm 3 RS_2 584
sql_rm 3 RS_2 585
sql_rm 3 RS_2 598
sql_rm 3 RS_2 599
sql_rm 3 RS_2 601
sql_rm 3 RS_2 602
sql_rm 3 RS_2 608
sql_rm 3 RS_2 616
sql_rm 3 RS_2 622
sql_rm 3 RS_2 624
sql_rm 3 RS_2 629
sql_rm 3 RS_2 632
sql_rm 3 RS_2 637
sql_rm 3 RS_2 651
sql_rm 3 RS_2 655
sql_rm 3 RS_2 661
sql_rm 3 RS_2 662
sql_rm 3 RS_2 667
sql_rm 3 RS_2 671
sql_rm 3 RS_2 681
sql_rm 3 RS_2 685
sql_rm 3 RS_2 694
sql_rm 3 RS_2 700
sql_rm 3 RS_2 711
sql_rm 3 RS_2 714
sql_rm 3 RS_2 728
sql_rm 3 RS_2 729
sql_rm 3 RS_2 734
sql_rm 3 RS_2 763
sql_rm 3 RS_2 766
sql_rm 3 RS_2 780
sql_rm 3 RS_2 784
sql_rm 3 RS_2 785
sql_rm 3 RS_2 789
sql_rm 3 RS_2 790
sql_rm 3 RS_2 797
sql_rm 3 RS_2 805
sql_rm 3 RS_2 829
sql_rm 3 RS_2 830
sql_rm 3 RS_2 839
sql_rm 3 RS_2 841
sql_rm 3 RS_2 844
sql_rm 3 RS_2 846
sql_rm 3 RS_2 856
sql_rm 3 RS_2 860
sql_rm 3 RS_2 863
sql_rm 3 RS_2 870
sql_rm 3 RS_2 872
sql_rm 3 RS_2 893
sql_rm 3 RS_2 896
sql_rm 3 RS_2 917
sql_rm 3 RS_2 931
sql_rm 3 RS_2 932
sql_rm 3 RS_2 934
sql_rm 3 RS_2 935
sql_rm 3 RS_2 938
sql_rm 3 RS_2 939
sql_rm 3 RS_2 941
sql_rm 3 RS_2 947
sql_rm 3 RS_2 960
sql_rm 3 RS_2 974
sql_rm 3 RS_2 977
sql_rm 3 RS_2 978
sql_rm 3 RS_2 982
sql_rm 3 RS_2 989
sql_rm 3 RS_2 992
