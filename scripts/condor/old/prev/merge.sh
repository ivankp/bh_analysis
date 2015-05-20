#! /bin/bash

BIN=/home/ivanp/work/bh_analysis/bin

mkdir -p ../hists/hists
mkdir -p ../hists/plots

all=""

for dir in H1j-ggf H2j-ggf H3j-ggf AA2j AA3j
do
  all+=" "`find ../$dir/hist -type d -name all`
done

for dir in $all; do
  #echo $dir;
  process=`echo $dir | cut -d'/' -f2`
  energy=`echo $dir | cut -d'/' -f4`

  namesuf=`echo $dir | cut -d'/' -f2,4- | rev | cut -d '/' -f 2- | rev | sed -e "s/\//\_/g"`;
  suffix="../hists/hists/$namesuf"
  echo $suffix

  for p in B RS V I; do
    hists="$suffix"_"$p".root
    if [ ! -f "$hists" ]; then
      if [ `ls $dir | grep '_'$p'_' | wc -l` -ne 0 ]; then
        hadd $hists $dir/*_"$p"_*.root
      elif [ `ls $dir | grep '^'$p'_' | wc -l` -ne 0 ]; then
        hadd $hists $dir/"$p"_*.root
      else
        echo "no $p hists for:" $dir
      fi
    fi
  done

  nlo="$suffix"_NLO.root
  if [ ! -f "$nlo" ]; then
    if [[ -f "$suffix"_B.root && -f "$suffix"_RS.root && -f "$suffix"_V.root && -f "$suffix"_I.root ]]; then
      $BIN/merge_parts $nlo "$suffix"_B.root "$suffix"_RS.root "$suffix"_V.root "$suffix"_I.root
    fi
  fi

  for p in B RS V I NLO; do
    hists="$suffix"_"$p".root
    plots="../hists/plots/$namesuf"_"$p".pdf
    if [[ ! -f "$plots" && -f "$hists" ]]; then
      if [[ $suffix == *"HT2-unc"* ]]; then
        $BIN/plot $hists -o $plots -t ": $process $p $energy"
      elif [[ $suffix == *"HT2-mH-cmp"* ]]; then
        $BIN/overlay $hists $plots
      fi
    fi
  done

done

