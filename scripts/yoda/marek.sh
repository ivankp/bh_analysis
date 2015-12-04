#!/bin/bash

path=`sed 's/marek\.sh$//' <<< $0`
fin=$1 # H3j_13TeV_NLO_HT2-unc_AntiKt4.yoda
fin_tmp=false

if [[ "$fin" == *.tar.bz2 ]]; then
  fin=/tmp/`tar xvjf $fin -C /tmp/`
  fin_tmp=true
fi

for dir in `sed -n 's|^Path=.*/\([^/:]*\)/[^/]*$|\1|p' $fin | uniq`; do
  fouta=`${path}marek.py $fin $dir`
  fout=`sed 's/_Jet[^_\.]*//;s/_PDF/_/' <<< $fouta`
  mv $fouta $fout
  echo $fout
  sed -i '
s|\([/=]\)H_pT|\1H_pT_incl|;
s|\([/=]\)H1j_pT|\1Hj_pT_incl|;
s|\([/=]\)H2j_pT|\1Hjj_pT_incl|;
s|\([/=]\)H3j_pT|\1Hjjj_pT_incl|;
s|\([/=]\)H4j_pT|\1Hjjjj_pT_incl|;
s|\([/=]\)jet\([0-9]\)_tau|\1tau_jet\2|;
s|\([/=]\)jets_tau_max|\1tau_jet_max|;
s|\([/=]\)jets_tau_sum|\1sum_tau_jet|;
s|\([/=]\)\(.*\)jjfb_|\1\2jjdy_|;
s|\([/=]\)jj_loose|\1loose|;
/_\(excl\|incl\)$/!s|\([/=]H_\?j\+.*\)|\1_incl|;
s/\(jet[0-9]_\(pT\|mass\)\)$/\1_incl/;
s/\([\/=]\)\([^\/=]*\)_d\(phi\|y\)/\1delta\3_\2/;
s/\([\/=]\)deltaphi_H_jj/\1deltaphi_Hjj/;
/[\/=]delta/s/\(H\?jj\)\(pT\)/\1/;
s/jets_HT/HT_jets/;
s/jets_N_\(excl\|incl\)/NJet_\1_30/;
s/H_jjdy_phi2_incl/deltaphi2/;
s/H_jj_phi2_incl/deltaphi2_jjpT/;
s/Hjjdy_mass_incl/H_dijet_mass_VBF/;
s/HjjpT_mass_incl/H_dijet_mass/;
s/deltaphi_jj$/deltaphi_jj_incl/;
s/deltaphi_jj_2j_excl/deltaphi_jj_excl/;
s/deltay_H_jj_incl/deltay_H_jj/;
s/HjjpT_mass_incl/H_dijet_mass/;
s/jjpT_mass/dijet_mass/
' $fout
done

if $fin_tmp; then
  rm $fin
fi

