#!/bin/bash

for f in `ls ~/disk2/hist/merged`; do
  name=~/disk2/hist/yoda/`basename $f .root`
  if [ ! -f "$name.yoda.tar.bz2" ]; then
    root2yoda ~/disk2/hist/merged/$f $name.yoda
    tar cjf $name.yoda.tar.bz2 $name.yoda
    rm $name.yoda
  fi
done
cd ..

