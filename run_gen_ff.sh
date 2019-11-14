#!/bin/sh

PARMCHK=parmchk2
# command -v $PARMCHK >/dev/null 2>&1 || { exit 1; }
if ! [ -x $(command -v $PARMCHK) ]
then
    PARMCHK=parmchk
    if ! [ -x $(command -v $PARMCHK) ]
    then
        echo -e "\n\tparmchk or parmchk2 not found. Quitting ..."
        exit 1
    fi
fi

rm -f ff.log
for i in abc; do 
   antechamber -i ${i}.pdb -fi pdb -o ${i}_GAFF.mol2 -fo mol2 -c bcc -nc 0 -at gaff2
   
   $PARMCHK -i ${i}_GAFF.mol2 -f mol2 -o ${i}.frcmod -a N

   j=$( Rscript -e "cat(toupper('$i'))" )
      
   cat >ff.leap << EOF
logFile leap.log
source leaprc.protein.ff14SB
source leaprc.gaff2
loadamberparams ${i}.frcmod
$j = loadmol2 ${i}_GAFF.mol2
saveoff $j ${i}.lib
quit
EOF
   
   tleap -s -f ff.leap &>> ff.log
done  

# To use the force field, 
# source leaprc.gaff2
#loadOff XX.lib 
#loadamberparams XX.frcmod
