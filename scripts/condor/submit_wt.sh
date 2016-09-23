#!/bin/bash

ana=/home/ivanp/work/bh_analysis
dir=/msu/data/t3work2/ivanp
out=$dir/condor

max_conc=150

mkdir -p $out

db="$dir/ntuples.db"

#####################################################################

for set in `sqlite3 $db "SELECT bh.id, wt.dir, wt.file FROM bh JOIN wt ON bh.id = wt.bh_id WHERE dset = 10 or dset = 11 or dset = 12"`; do

while [ "`condor_q ivanp | tail -1 | sed 's/\([0-9]*\) jobs.*/\1/'`" -ge "$max_conc" ]; do sleep 1800; done

arr=(`echo $set | tr '|' ' '`)
wt="${arr[1]}/${arr[2]}"

mkdir -p ${arr[1]}

if [ -a "$wt" ]; then
  echo "."
  continue;
fi

name=`basename ${arr[2]} .root`
echo $name

bh=`sqlite3 -separator '/' $db "SELECT dir, file FROM bh WHERE id = ${arr[0]}"`

cmd="$ana/bin/reweigh --bh=\"$bh\" -c \"$ana/config/rew_Ht2.xml\" -o \"$wt\" --pdf=\"CT10nlo\" --tree-name=\"t3\""

# Form temporary wrapper script

echo "#!/bin/bash" > $out/$name.sh
#echo "LD_LIBRARY_PATH=\`echo \$LD_LIBRARY_PATH | sed 's|\(/lib64:/usr/lib64\):\(.*\)|\2:\1|'\`" >> $out/$name.sh
echo 'LD_LIBRARY_PATH=/home/ivanp/local/gcc/lib64:/home/ivanp/local/lib:$LD_LIBRARY_PATH' >> $out/$name.sh
echo 'echo $LD_LIBRARY_PATH' >> $out/$name.sh
echo "ldd $ana/bin/reweigh" >> $out/$name.sh
echo "$cmd" >> $out/$name.sh
chmod +x $out/$name.sh

# Form condor script

echo "
Universe   = vanilla
Executable = $out/$name.sh
Error      = $out/$name.err
Output     = $out/$name.out
Log        = $out/$name.log
getenv = True
Queue
" | condor_submit - > /dev/null

sleep 1

done

echo 'DONE!'

