#!/bin/bash
## Have all files in the same folder that you are running the script in
## Set all values below to your specification then run script by using source new_mmpbsa.sh
## No need to put .prmtop after file

prmtop_complex='sys_lig_nowat'
prmtop_receptor='sys_apo_nowat'
prmtop_ligand='lig_nowat'
trajin_complex='complex_mmpbsa.nc'
trajin_receptor='receptor_mmpbsa.nc'
trajin_ligand='ligand_mmpbsa.nc'
# Colon must be placed before entering ranges to be stripped
#receptor_strip=':165-999999'
receptor_strip=':9999'
#ligand_strip=':1-164'
ligand_strip=':9999'

mkdir complex
mkdir receptor
mkdir ligand
mkdir rstfiles


initialize_files () {
    fin=trajin_$1
    cat > collect_$1.ptraj <<EOF
trajin ${!fin}
$4 $5
trajout rstfiles/$1 restart multi
go
EOF


#    fprmtop=prmtop_$1.prmtop
    cat > $2.sh <<EOF
#!/bin/bash 
export LD_LIBRARY_PATH="/software/usr/cuda-8.0/lib64:/usr/lib64/openmpi/lib/:\${LD_LIBRARY_PATH}"
export CUDA_VISIBLE_DEVICES="0"

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

cpptraj $3.prmtop collect_$1.ptraj
EOF

    cat >$1.sh<<EOF
#!/bin/bash 
export LD_LIBRARY_PATH="/software/usr/cuda-8.0/lib64:/usr/lib64/openmpi/lib/:\${LD_LIBRARY_PATH}"
export CUDA_VISIBLE_DEVICES="0"

echo Time is `date`
echo Directory is `pwd`
echo Reading from \$1 to \$2

date
x=\$1  
total=\$2

while [ \$x -le \$total ]
do

    sander -O -i ../../mmpbsa$6.in -o mmpbsa_c.out -p ../../$3.prmtop -c ../../rstfiles/$1.\${x}
    awk 'NR==6, NR==9 {print}' mdinfo | sed -e 's/1-4 /1-4-/g' | awk '{print \$3,  \$6,  \$9}' | xargs   >> energies_$1.dat

    x=\$((\$x + 1 ))
done
date
EOF

    chmod +x $1.sh
    chmod +x $2.sh

    nohup sh $2.sh &
}

initialize_files complex c ${prmtop_complex}
initialize_files receptor r ${prmtop_receptor} strip ${receptor_strip}
initialize_files ligand l ${prmtop_ligand} strip ${ligand_strip} _l

cat >mmpbsa.in <<EOF
initial minimization w/ position restraints on DNA, 9.0 cut
 &cntrl
  imin   = 1, ntx  = 1, ipb  = 2, inp = 2, 
 /
 &pb
  epsout=80.0, epsin=1.0, fillratio=2.0,
  radiopt = 0,
 /
EOF

cat >mmpbsa_l.in <<EOF
initial minimization w/ position restraints on DNA, 9.0 cut
 &cntrl
  imin   = 1, ntx  = 1, ipb  = 2, inp = 2, 
 /
 &pb
  epsout=80.0, epsin=1.0, fillratio=4.0,
  radiopt = 0,
 /
EOF
