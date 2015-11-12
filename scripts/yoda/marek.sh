#!/bin/bash

fin=H3j_13TeV_NLO_HT2-unc_AntiKt4.yoda

for dir in `sed -n 's|^Path=.*/\([^/:]*\)/[^/]*$|\1|p' $fin | uniq`; do
  fouta=`./scripts/marek.py $fin $dir`
  fout=`sed 's/_Jet[^_\.]*//;s/_PDF/_/' <<< $fouta`
  mv $fouta $fout
  echo $fout
  sed -i '
s|\([/=]\)H_pT|\1H_pT_incl|;
s|\([/=]\)H1j_pT|\1H_j_pT_incl|;
s|\([/=]\)H2j_pT|\1H_jj_pT_incl|;
s|\([/=]\)H3j_pT|\1H_jjj_pT_incl|;
s|\([/=]\)H4j_pT|\1H_jjjj_pT_incl|;
s|\([/=]\)jet\([0-9]\)_tau|\1tau_jet\2|;
s|\([/=]\)jets_tau_max|\1tau_jet_max|;
s|\([/=]\)jets_tau_sum|\1sum_tau_jet|;
s|\([/=]\)jets_tau_sum|\1sum_tau_jet|
' $fout
done


