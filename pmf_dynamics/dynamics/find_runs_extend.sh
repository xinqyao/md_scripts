#!/bin/bash

root=ts_time_life1
maxrun=100

for i in $root/*; do
   j=$(basename $i .txt | sed 's/.*_//')
   ofile=to_extend_$j
   Rscript -e "dat <- read.table('$i')" \
           -e "dat <- dat[dat[, 3]!=-1, ]" \
           -e "cat(setdiff(seq(1, $maxrun), as.numeric(dat[, 1])), sep='\n')" \
           > $ofile
done

## rerun
for i in $root/*; do
   j=$(basename $i .txt | sed 's/.*_//')
   ofile=to_rerun_$j
   Rscript -e "dat <- read.table('$i')" \
           -e "dat <- dat[dat[, 3]!=0, ]" \
           -e "cat(setdiff(seq(1, $maxrun), as.numeric(dat[, 1])), sep='\n')" \
           > $ofile
done
