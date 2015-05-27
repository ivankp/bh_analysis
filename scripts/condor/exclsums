#!/bin/sh

dir=~/disk2/hist/merged
BIN=../../bin

for p in "H" "AA"; do
for n in 2 3; do
for E in 8 13; do

for f in `ls $dir | grep "$p$n"j_"$E"TeV | grep _NLO_`; do

  props=`echo $f | sed 's/^[^_]*_[^_]*_[^_]*_\(.*\)/\1/'`
  #echo $props

  bad=0
  args=()
  for t in "RS" "V" "I" "NLO"; do

    nj=$n
    if [ "$t" == "NLO" ]; then ((--nj)); fi

    arr=(`ls $dir | grep "$p$nj"j_"$E"TeV_"$t[^_]*"_"$props"`)
    if [ "${#arr[@]}" -eq "0" ]; then
      #echo "bad t is $t"
      ((++bad))
      break
    else
      for ((i=0; i < ${#arr[@]}; i++)); do
        arr[$i]="$dir/${arr[$i]}"
      done
      args+=(${arr[@]})
    fi

  done

  if [ "$bad" -gt "0" ]; then continue; fi

  es=$dir/$p$((n-1))"j_"$E"TeV_ES_"$props

  if [ ! -f "$es" ]; then
    $BIN/merge_parts $es ${args[@]}
  fi

done

done
done
done

