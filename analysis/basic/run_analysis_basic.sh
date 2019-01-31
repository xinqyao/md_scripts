#!/bin/bash

## !!!! This script is supposed to run under $TRAJ_ROOT/analysis/basic !!!

#imin=0
#if test $# -ge 1 && test $1 == "em"; then imin=1; fi

perl process_mdout.perl `ls ../../*.out -t -r`

## skip all lines with #col not equal to 2
for i in summary.*; do
   awk 'NF==2{print}' $i > tmp
   mv tmp $i
done

#if [ $imin -ge 1 ]; then
#   mv summary.VOLUME summary.VOLUME.ori
#   mv summary.DENSITY summary.DENSITY.ori
#   sed '1,11'd summary.VOLUME.ori > summary.VOLUME
#   sed '1,11'd summary.DENSITY.ori > summary.DENSITY
#fi

#ptraj $WP/md/sys_gdp_nowat.prmtop <<EOF
#trajin $WPmd/dyna.01.mdcrd
#reference $WP/md/sys_gdp_nowat.inpcrd 
##center :1-314 mass origin
##image origin center
##strip :WAT
#rms reference out bb.rms :1-314@N,CA,C,O time 2.0
#rms reference out ca.rms :1-314@CA time 2.0
#rms reference out all.rms :1-314 time 2.0
#EOF
# for R
Rscript plot_analysis_basic.r
