#!/bin/bash

## only used for peptides in complex with an enzyme

if test $# -lt 1; then
   echo Usage: $0 para
   exit 0
fi

## Full length #########
#atom1=":63@CE"
#atom2=":167@P"

## cutoff 10 angstrom
atom1=":157@CG"
atom2=":168@CG"
nres=170
########################

## PPIASE DOMAIN ONLY ##
#atom1=":13@CE"
#atom2=":117@P"

## cutoff 10 angstrom
#atom1=":107@CG"
#atom2=":118@CG"
#nres=113
#####################E##


para=$1
shift
root=$(pwd)
top=($(ls $root/$para/*.prmtop))
#top=$(basename $top)

if ! test -r check_binding/$para; then
   mkdir -p check_binding/$para
fi

cd check_binding/$para
files=( $(ls $root/$para/production.-[0-9][0-9][0-9].restrt 2>/dev/null | sort -rn) \
        $(ls $root/$para/production.-[0-9][0-9].restrt 2>/dev/null | sort -rn) \
        $(ls $root/$para/production.-[0-9].restrt 2>/dev/null | sort -rn) \
        $(ls $root/$para/production.[0-9].restrt 2>/dev/null | sort -n) \
        $(ls $root/$para/production.[0-9][0-9].restrt 2>/dev/null | sort -n) \
        $(ls $root/$para/production.[0-9][0-9][0-9].restrt 2>/dev/null | sort -n)\
)
angs=($(echo ${files[@]} | xargs -n 1 basename | sed 's/production.\([-0-9]*\).restrt/\1/'))
echo parm $top > ptraj.in
for((i=0; i<${#files[@]}; i++)); do
#for((i=${#files[@]}-1; i>=0; i--)); do
   echo trajin ${files[i]} >> ptraj.in
done
echo "strip :WAT,Na+,Cl-" >> ptraj.in
echo "center :1-$nres" >> ptraj.in
echo image center familiar >> ptraj.in
echo "rms first :1-$nres@N,CA,C,O" >> ptraj.in
echo "trajout restrt.pdb" >> ptraj.in
echo "distance out dist.out $atom1 $atom2" >> ptraj.in
echo go >> ptraj.in
cpptraj -i ptraj.in &> err.log


echo "dat <- read.table('dist.out')" > plot.r
echo "png(file='plot.png')" >> plot.r
echo "plot(x=as.numeric(strsplit('${angs[@]}', split='\\\\s+')[[1]])," >> plot.r
echo "     y=dat[, 2], xlab='Angle', ylab='Distance ($atom1 - $atom2; A)', typ='o')" >> plot.r
echo "abline(h=10, lty=3, col='red', lwd=2)" >> plot.r
echo "dev.off()" >> plot.r
Rscript plot.r

cd $root
echo Done.

