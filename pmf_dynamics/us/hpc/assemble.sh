#!/bin/bash

para='14.0'
root='../..'
bd=FALSE

################
### obsolete ###
## Randomly choose 5 for checking. It will be common for both us and dynamics
#tmploc=$(find . -name production.5.restrt)
#ref=$root/${para}/production.5.restrt
###############

## Choose the first available .restrt for checking. It will be common for both us and dynamics
tmploc=( $(find . -name production.[-0-9]*.restrt  | grep -e ^"./${para}_" -e ^"${para}_") )
ref=$root/${para}/$(basename $tmploc)

if test -r $ref && diff $tmploc $ref &> /dev/null; then
   echo Seems you have already done the assembly...stopped
   exit 0
fi


echo Copy files to $root/${para}
echo -n 'Continue? (y/n)'
read a
if test ! $a == "y" && test ! $a == "Y"; then
   exit 0
fi

if test $bd == "TRUE"; then
   if ! test -d $root/$para/restrt; then
      rsync -av $para/ $root/$para/
   fi
   exit 0
fi

echo Have you trimmed? '(Always "yes" if not initial replica dynamics or extending us or dynamics)'
echo -n 'Continue? (y/n)'
read b
if test ! $b == "y" && test ! $b == "Y"; then
   exit 0
fi

if test $a == "y" || test $a == "Y"; then
   for i in ${para}_*; do
      angs=( $i/ang.[-0-9]*.out ) 
      for j in ${angs[@]}; do
        if test -r ${root}/${para}/$(basename $j); then
           cat $j >> ${root}/${para}/$(basename $j)
        else 
           rsync -av $j ${root}/${para}/
        fi
      done
      if test -r ${root}/${para}/restrt; then
         rsync -av $i/ --exclude 'ang.[-0-9]*.out' --exclude 'restrt' ${root}/${para}/
      else
         rsync -av $i/ --exclude 'ang.[-0-9]*.out' ${root}/${para}/
      fi
   done
fi
