#!/bin/bash

rm *dat
rm *mmpbsa

if [ "$#" -ne 1 ]; then
    echo "Enter output file name"
    exit
fi

i=0
j=0
k=0
#`
for i in $(seq 10);do
cat complex/complex$i/*.dat >>energies_complex.dat
done

for j in $(seq 10);do
    cat ligand/ligand$j/*.dat >>energies_ligand.dat
    ((j+=1))
done

for k in $(seq 10);do
   cat receptor/receptor$k/*.dat >>energies_receptor.dat
   ((k+=1))
done

awk '{print ($1 + $2 + $3 + $5 + $8 + $4 + $7 + $6 + $10 + $11)}'  energies_complex.dat > PBTOT_complex.dat
awk '{print ($1 + $2 + $3 + $5 + $8 + $4 + $7 + $6 + $10 + $11)}'  energies_receptor.dat > PBTOT_receptor.dat
awk '{print ($1 + $2 + $3 + $5 + $8 + $4 + $7 + $6 + $10 + $11)}'  energies_ligand.dat > PBTOT_ligand.dat

#paste PBTOT_complex.dat PBTOT_receptor.dat PBTOT_ligand.dat > energies_c_r_l.dat
#awk '{ print $1 - $2 - $3}' energies_c_r_l.dat > MMPBSA_$1.mmpbsa

Rscript final.r
