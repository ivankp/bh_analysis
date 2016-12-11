#!/bin/bash

# ./bin/hist_Hjets_mtop_study -j2 -s config/Hjets_basic.css -o inf_raw.root \
#   --bh=/msu/data/t3work4/luisonig/H2jets_ggf/NTuplesFiles/H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_100.root
# ./bin/merge_parts inf.root inf_raw.root

./bin/hist_Hjets_mtop_study -j2 -s config/Hjets_basic.css -o finite_raw.root \
  --has-ncount \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_101.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_102.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_103.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_104.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_105.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_106.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_107.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_108.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_109.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_110.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_111.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_112.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_114.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_116.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_117.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_119.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_120.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_121.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_122.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_123.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_124.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_125.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_126.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_127.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_128.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_129.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_130.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_131.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_132.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_133.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_134.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_135.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_136.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_137.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_138.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_139.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_140.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_141.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_142.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_143.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_144.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_145.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_146.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_147.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_148.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_149.root \
  --bh=/msu/data/t3work4/luisonig/H2jets_ggf_mtop/NTuplesFiles/RWGT1_H2.0j_GGFHT_B_6500_pt25.0_eta4.5_r100_150.root
./bin/merge_parts finite.root finite_raw.root

