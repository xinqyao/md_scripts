#!/bin/bash

#for i in 4.9 5.6; do
#for i in 14.0; do
#for i in 4.9 5.6 6.3 14.0; do
for i in 3.5 4.2; do
#for i in 4.9; do
   echo Processing ${i}...
   cd $i
#   ./run.sh 321
   ./run.sh 181
#   ./run.sh 71
#   ./run.sh
   cd ..
done

echo Done.
