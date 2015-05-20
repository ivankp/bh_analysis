#!/bin/sh

BIN=/home/ivanp/work/bh_analysis/bin

for E in 8 13; do
  $BIN/merge_parts ../hists/exclsums/H2-3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.root \
    ../hists/hists/H3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_RS.root \
    ../hists/hists/H3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_I.root \
    ../hists/hists/H3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_V.root \
    ../hists/hists/H2j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_B.root \
    ../hists/hists/H2j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_RS.root \
    ../hists/hists/H2j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_I.root \
    ../hists/hists/H2j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4_V.root

  $BIN/plot -i ../hists/exclsums/H2-3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.root \
            -o ../hists/exclsums/H2-3j-ggf_"$E"TeV_HT2-unc_AntiKt4_CT10nlo_pt30-eta4.4.pdf \
            -t ": H2&3j ExclSum "$E"TeV"
done

