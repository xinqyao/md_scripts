#!/bin/bash

if test $# -lt 2; then
   echo Usage: $0 para N_ns
   exit 0
fi

para=$1
shift
root=$(pwd)

# how many chunks?
#nchunk=120
nchunk=$1
shift
nline=$(($nchunk * 50000))
#for i in 3.5 4.2 4.9 5.6; do
#for i in 14.0; do
#for i in $para; do
#   wc -l $i/ang.[-0-9]*.out | awk "!/total/&&\$1!=${nline}{print $i, \$0}"
#done
for((i=-175; i<=180; i+=5)); do
   if [[ $i -gt -150 ]] && [[ $i -lt -30 ]]; then   # sampling range -30 ~ 210 (-150)
      continue
   fi
   file=$para/ang.${i}.out
   if ! test -r $file; then
      echo $file missing 
   else
      wc -l $para/ang.${i}.out | awk "!/total/&&\$1!=${nline}{print $i, \$0}"
   fi
done
