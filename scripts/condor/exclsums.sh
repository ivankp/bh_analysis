#!/bin/sh

BIN=../../bin

for n in 1 2; do
m=$(($n + 1))
for E in 8 13; do
  hadd exclsums/H"$n$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.root \
    hists/H"$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_RS.root \
    hists/H"$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_I.root \
    hists/H"$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_V.root \
    hists/H"$n"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_NLO.root

  $BIN/plot -i exclsums/H"$n$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.root \
            -o exclsums/H"$n$m"j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.pdf \
            -t ": H"$n"&"$m"j ExclSum "$E"TeV"
done
done
