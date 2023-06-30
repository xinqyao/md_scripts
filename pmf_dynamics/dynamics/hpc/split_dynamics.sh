#!/bin/bash

## Brownian dynamics?
bd=FALSE

para=6.3
root=$(pwd)/../..
prmtop_in=ffpspr_ppiase
prmtop=ffpspr

if test $bd == "TRUE"; then
   rsync -av $root/$para ./
   ii=01
   sed -e "s/%ii%/${ii}/g" -e "s/%dir%/${para}/g"  \
       -e "s/%prmtop_in%/${prmtop_in}/g"\
       -e "s/%prmtop%/${prmtop}/g"\
     rundynamics_bd_acids_template.sh > job${ii}.sh
   chmod +x job*.sh
   exit 0
fi

## Number of jobs; ignored if bd=TRUE
n=24
#n=4
#n=4

## Number of replicas in each job; ignored if bd=TRUE
k=2
#k=3
#k=2

## if bExtend = FALSE
ibase=0  # what is the last previous job ID?
abase=0 # what is the last previous "number"?
nmax=100 # the maximum number of replicas   


## Extended simulations?
bExtend=TRUE

## Rerun simulations?
bRerun=FALSE

## if bExtend = TRUE
feruns=$root/to_extend_$para # filename contain run No.'s to extend

## if bRerun = TRUE
frruns=$root/to_rerun_$para # filename contain run No.'s to rerun

if test $bExtend == TRUE || test $bRerun == TRUE; then

   if test $bExtend == TRUE; then
      eruns=( $(cat $feruns) )
   else
      eruns=( $(cat $frruns) )
   fi
#   echo ${eruns[@]}

   if test ${#eruns[@]} -eq 0; then
      echo No runs to extend or rerun.
      exit 0
   fi

#   ## clear up root folder
#   if test $bRerun == TRUE; then
#      for((i=0; i<${#eruns[@]}; i++)); do
#         rm $root/$para/ang.${eruns[i]}.out
#      done
#   fi

   kk=0
   bStop=FALSE
   for((j=1; j<=$k; j++)); do
      for((i=1; i<=$n; i++)); do
         if test $bStop == TRUE; then
            break 2 
         fi
         ii=$(echo $i | awk '{printf "%02d", $1}')
         if ! test -r ${para}_$ii; then
            mkdir -p ${para}_$ii/restrt
            cp $root/$para/*.prmtop ${para}_$ii/
         fi
         if test $bExtend == TRUE; then
            cp $root/$para/production.${eruns[kk]}.restrt ${para}_$ii/restrt/production.bd.${eruns[kk]}.restrt
         else
            ## clear up root folder
            rm $root/$para/ang.${eruns[kk]}.out
            cp $root/$para/restrt/production.bd.${eruns[kk]}.restrt ${para}_$ii/restrt/production.bd.${eruns[kk]}.restrt
         fi
         let kk=kk+1
         if test $kk -eq ${#eruns[@]}; then
            bStop=TRUE
         fi
      done
   done

   for i in $(ls -d ${para}_*/); do
      ii=$(basename $i | sed 's/.*_//')
      runs=( $(ls ${para}_$ii/restrt | sed 's/.*bd.\([0-9]*\).restrt/\1/') )
      a1=0
      a2=$(expr ${#runs[@]} - 1)
      sed -e "s/%ii%/$ii/g" -e "s/%dir%/${para}_$ii/g"  \
          -e "s/%a1%/$a1/g" -e "s/%a2%/$a2/g" \
          -e "s/%numbers%/($(echo ${runs[@]}))/g" \
        rundynamics_acids_template.sh > job${ii}.sh
      chmod +x job${ii}.sh
   done

   if test $kk -lt ${#eruns[@]}; then
      echo Done but not all runs are deployed.
   else
      echo Done with all runs deployed.
   fi
   
else

   for((i=1; i<=$n; i++)); do
      ii=$(echo $(expr $i + $ibase) | awk '{printf "%02d", $1}')
   #   echo $ii
      if ! test -r ${para}_$ii; then
         mkdir -p ${para}_$ii/restrt
      fi
      cp $root/$para/*.prmtop ${para}_$ii/
#      cp $root/$para/*.restrt ${para}_$ii/
      a1=$( echo $i $k $abase | awk '{print ($1-1)*$2+1+$3}' )
      a2=$( expr $a1 + $k - 1 )
      if test $a2 -gt $nmax; then
         a2=$nmax
      fi
      for((j=$a1; j<=$a2; j++)); do
         cp $root/$para/restrt/production.bd.${j}.restrt ${para}_$ii/restrt/
      done 
      sed -e "s/%ii%/$ii/g" -e "s/%dir%/${para}_$ii/g"  \
          -e "s/%a1%/$a1/g" -e "s/%a2%/$a2/g" \
          -e "s/%numbers%/(\$(seq 0 $a2))/g" \
       rundynamics_acids_template.sh > job${ii}.sh
      chmod +x job${ii}.sh
   done

fi
