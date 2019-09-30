#module load research
#module load pdb2pqr
#module load pymol/1.8.6

#mkdir prep
mkdir equil
#mkdir windows
#mkdir us
#mkdir dynamics

cd prep

cp /software/md_scripts_hamelberglab/frcmod.ff99SB_w_dih .

# use leap to build the peptide structures: 
# 1. Ace-Phe-Phe-Ser-Pro-Arg-Nme
# 2. Ace-Phe-Phe-Sep-Pro-Arg-Nme
#
# NOTE:
# - assume the -2 net charge form of phosphserine (SER -> SEP)
# - use our modified force field files
cp /software/md_scripts_hamelberglab/frcmod.phosaa10.ff14sb .
cp /software/md_scripts_hamelberglab/leaprc.phosaa10.ff14sb .
cp /software/md_scripts_hamelberglab/phos_amino94_ff14sb.lib .

# modify the parameters for the proline in the peptide
# compare to $AMBERHOME/dat/leap/parm/parm10.dat
# frcmod.mod

cat >tleap.in << EOF
source leaprc.protein.ff14SB
#source leaprc.phosaa10.ff14sb
source leaprc.water.tip3p
addAtomTypes{{"DN" "N" "sp2"}}
loadamberparams frcmod.ff99SB_w_dih
loadamberparams frcmod.mod
set default PBradii mbondi2

mol = sequence {ACE PHE PHE SER PRO ARG NME}

set mol.5.N type "DN"
savepdb mol ffspr.pdb

solvatebox mol TIP3PBOX 10.0 iso
addions mol Na+ 0
addions mol Cl- 0
savepdb mol ffspr_box.pdb
saveamberparm mol ffspr_box.prmtop ffspr_box.inpcrd

quit
EOF
tleap -f tleap.in &> log

cd ..

#cp /software/md_scripts_hamelberglab/new_equil.in .

## modify the script and execute it
vi new_equil.in
./new_equil.in

cd windows
cp /software/md_scripts_hamelberglab/umbrella/*.in ./
sed 's/8.33250000E+00/3.50000000E+00/' ../prep/ffspr_box.prmtop > ./ffspr_box.prmtop
## modify scripts
vi equil.in
#...
./umbrella_startup.in

cd ..

cd us
## create a file: runus.sh
## also, metafile
vi runus.sh

cd ..

#transfer files and execute one by one equil, windows, and us
# cd equil; ./jobe.in &> err.log
# cd ../windows; ./run_umbrella.in
# cd ../us; ./runus.sh
