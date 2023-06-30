#!/bin/bash

n=24
para=14.0
root=$(pwd)/../..

#files=( $(ls -rt $root/$para/production.[-0-9]*.restrt) )

## -175 to 180
#files=( $(ls $root/$para/production.[0-9][0-9][0-9].restrt | sort -rn) \
#        $(ls $root/$para/production.[0-9][0-9].restrt | sort -rn) \
#        $(ls $root/$para/production.[0-9].restrt | sort -rn) \
#        $(ls $root/$para/production.-[0-9].restrt | sort -n) \
#        $(ls $root/$para/production.-[0-9][0-9].restrt | sort -n) \
#        $(ls $root/$para/production.-[0-9][0-9][0-9].restrt | sort -n)\
#)

## -30 to 180, -175 to -150
files=( $(ls $root/$para/production.[0-9][0-9][0-9].restrt | sort -rn) \
        $(ls $root/$para/production.[0-9][0-9].restrt | sort -rn) \
        $(ls $root/$para/production.[0-9].restrt | sort -rn) \
        $(ls $root/$para/production.-[0-9].restrt | sort -n) \
        $(ls $root/$para/production.-[0-2][0-9].restrt | sort -n) \
        $(ls $root/$para/production.-30.restrt | sort -n) \
        $(ls $root/$para/production.-1[5-9][0-9].restrt | sort -n)\
)
#k=$(expr ${#files[@]} / $n )
while [ ${#files[@]} -gt 0 ]; do
   for((i=1; i<=$n; i++)); do
      ii=$(echo $i | awk '{printf "%02d", $1}')
   #   echo $ii
      if ! test -r ${para}_$ii; then
         mkdir ${para}_$ii
         cp $root/$para/*.prmtop ${para}_$ii/
      fi
      cp ${files[0]} ${para}_$ii/
      files=( ${files[@]:1} )
      if [ ${#files[@]} -eq 0 ]; then
         break
      fi
   done
done
for((i=1; i<=$n; i++)); do
   ii=$(echo $i | awk '{printf "%02d", $1}')
   angles=( $(ls ${para}_$ii/*.restrt | xargs -n 1 basename | sed 's/production.\([-0-9]\+\).restrt/\1/') )
   a1=0
   a2=$(expr ${#angles[@]} - 1)
   sed -e "s/%ii%/$ii/g" -e "s/%dir%/${para}_$ii/g"  \
       -e "s/%a1%/$a1/g" -e "s/%a2%/$a2/g" \
       -e "s/%angles%/($(echo ${angles[@]}))/g" \
       runus_acids_template.sh > job${ii}.sh
   chmod +x job${ii}.sh
done
