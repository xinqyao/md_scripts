#!/bin/bash

# how many chunks?
nchunk=50
nline=$(($nchunk * 50000))
for i in 3.5 4.2 4.9 5.6; do
   wc -l $i/ang.[-0-9]*.out | awk "!/total/&&\$1!=${nline}{print $i, \$0}"
done
