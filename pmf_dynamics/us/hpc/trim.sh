#!/bin/bash

if test $# -lt 1; then 
   echo WARNING: this script will cut out the first line of data!!
   echo Usage: ./trim.sh go
   exit 1
fi

## somehow amber output the time 0 data twice if istrt=0
## remove the redundant top line
for i in [0-9]*_*/ang.[0-9]*.out; do
   echo $i
   sed -i "1"d $i
done
echo Done.
