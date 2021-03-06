#!/bin/bash

dir=/home/ivanp/disk2/hist

mkdir -p $dir/merged
mkdir -p $dir/plots

cnt1=0
cnt2=0
cnt3=0

#if false; then

# hadd same parts -----------------------------------------
for f in `ls "$dir/raw" |
  sed 's/\(.*\)_[0-9]\+-[0-9]\+\.root$/\1/' | sort | uniq`
do

  arr=(`ls "$dir/raw" | grep "^${f}_[0-9]*-[0-9]*\.root$"`)
  for ((i=0; i < ${#arr[@]}; i++)); do
     arr[$i]=`echo "$dir/raw/${arr[$i]}"`
  done

  if [ ! -f "$dir/merged/$f.root" ]; then
    if [ "${#arr[@]}" -eq "1" ]; then
      # copy if there's only one file
      cp -v ${arr[0]} "$dir/merged/$f.root"
    else
      # merge if there are many
      hadd "$dir/merged/$f.root" ${arr[*]}
    fi
    ((++cnt1))
  fi

done
if [ "$cnt1" -eq "0" ]; then echo "Nothing to hadd"; fi

#fi

# merge different parts -----------------------------------
for f in `ls "$dir/merged" | grep -v "_NLO_" | grep -v "_ES_" |
  sed -n 's/\(.*\)_\(\(B\|I\|RS\|V\)\(-[^_]*\)\?\)_\(.*\)/\1\.\*\5/p' | sort | uniq`
do

  arr=(`ls "$dir/merged" | grep "^$f" | grep -v "_NLO_"`)
  for ((i=0; i < ${#arr[@]}; i++)); do
     arr[$i]=`echo "$dir/merged/${arr[$i]}"`
  done

  nlo="$dir/merged/"`echo $f | sed 's/\.\*/_NLO_/'`
  if [ ! -f "$nlo" ]; then
    ../../bin/merge_parts $nlo ${arr[*]}
    echo ""
    ((++cnt2))
  fi

done
if [ "$cnt2" -eq "0" ]; then echo "Nothing to merge"; fi

#exit

# make plots ----------------------------------------------
for f in `ls "$dir/merged"`; do

  plot="$dir/plots/`basename $f .root`.pdf"

  if [ ! -f "$plot" ]; then
    if [[ $f == *"HT2-mH-cmp"* ]]; then
      ../../bin/overlay "$dir/merged/$f" "$plot"
    else
      exe="../../bin/plot $dir/merged/$f -o $plot"
      exe=$exe`echo $f | sed 's/\([^_]*\)_\([^_]*\)_\([^_]*\).*/ -t \": \1 \3 \2\"/'`
      exe=$exe`echo $f | sed -n 's/.*VBF\([^_\.]*\).*/ -l \"VBF cut: \1\"/p'`
      exe=$exe`echo $f | sed -n 's/.*MAA\([^_\.]*\).*/ -l \"m_\{AA\} cut: \1\"/p'`
      eval $exe
    fi
    echo ""
    ((++cnt3))
  fi

done
if [ "$cnt3" -eq "0" ]; then echo "Nothing to plot"; fi

echo "DONE!"

