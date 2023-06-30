#!/bin/bash
module load wham/2.0.9.1

sep=1
para=3.5

maxn=$(cat ../../../$para/ang.0.out | wc -l | awk '{print $1/50000}')
if echo $maxn | grep -q "\\."; then 
   echo Number of cycles not an integer.
   exit 1
fi

if test $# -eq 0; then
   bb=1
   ee=$maxn
elif test $# -eq 1; then
   bb=$1
   ee=$maxn
else
   bb=$1
   ee=$2
   if test $ee -gt $maxn; then 
      echo Number of cycles greater than the maximum.
      exit 1
   fi
fi

rm -rf tmp
mkdir tmp

for ((rnum=$bb; rnum<=$ee; rnum++)); do
      for ifile in ../../../$para/ang.{-30,-[12][0-9],-[0-9],[0-9]*}.out; do
         echo $para \(sep=$sep\), $rnum \($maxn\), $(basename $ifile .out)
         if test $bb -gt 1 && test $rnum -eq $bb; then
            awk "(NR>0) && (NR<=50000*($bb-1)) && ((NR-1) % $sep==0){a=\$2; print \$1, a}" $ifile > tmp/$(basename $ifile)
         fi
         awk "(NR>50000*($rnum-1)) && (NR<=50000*$rnum) && ((NR-1) % $sep==0){a=\$2; print \$1, a}" $ifile >> tmp/$(basename $ifile)
      done
      for ifile in ../../../$para/ang.-1[5-7][0-9].out; do
         echo $para \(sep=$sep\), $rnum \($maxn\), $(basename $ifile .out)
         if test $bb -gt 1 && test $rnum -eq $bb; then
            awk "(NR>0) && (NR<=50000*($bb-1)) && ((NR-1) % $sep==0){a=\$2; a=a+360; print \$1, a}" $ifile > tmp/$(basename $ifile)
         fi
         awk "(NR>50000*($rnum-1)) && (NR<=50000*$rnum) && ((NR-1) % $sep==0){a=\$2; a=a+360; print \$1, a}" $ifile >> tmp/$(basename $ifile)
      done

#      wham P -63 297 100 0.000001 300.0 0 metafile pmf.out 10 1234 &> err.log
      wham -31.5 208.5 80 0.000001 300.0 0 metafile pmf.out 10 1234 &> err.log
      awk 'NR > 1 && NR <82{print $1,$2,$3}' pmf.out > fep.${rnum}.out
done
rm -rf tmp
