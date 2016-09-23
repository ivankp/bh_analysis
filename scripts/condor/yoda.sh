#!/bin/bash

cd ~/disk2/hist/yoda/
for f in `ls ../merged`; do
  name=`basename $f .root`
  if [ ! -f "$name.yoda.tar.bz2" ]; then
    root2yoda ../merged/$f $name.yoda
    tar cjf $name.yoda.tar.bz2 $name.yoda
    rm $name.yoda
  fi
done
cd -

