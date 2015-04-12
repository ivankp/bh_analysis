#! /bin/bash

cat out/*.err

echo `ls out/*.err | wc -l` files

if [ "$1" == "hist" ]; then
  all_ok=true
  for f in out/*.out; do
    if [[ "`tail -1 $f`" != "Processed events"* ]]; then
      echo "$f is incomplete"
      tail -1 $f
      echo ""
      all_ok=false
    fi
  done
  if [ "$all_ok" == true ]; then
    echo "All jobs complete"
  fi

elif [ "$1" == "wt" ]; then
  all_ok=true
  for f in out/*.out; do
    if [[ "`tail -1 $f`" != *"Wrote"* ]]; then
      echo "$f is incomplete"
      tail -1 $f
      echo ""
      all_ok=false
    fi
  done
  if [ "$all_ok" == true ]; then
    echo "All jobs complete"
  fi

fi

