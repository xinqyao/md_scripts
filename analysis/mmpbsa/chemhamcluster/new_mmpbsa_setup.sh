#!/bin/bash
## Have all files in the same folder that you are running the script in
## Set all values below to your specification then run script by using source new_mmpbsa.sh
## No need to put .prmtop after file

prmtop_complex='pcbw_nowt'
prmtop_receptor='pcbw_receptor'
prmtop_ligand='b'
trajin='pcbw.nc'
# Colon must be placed before entering ranges to be stripped
receptor_strip=':158-999999'
ligand_strip=':1-157'

mkdir complex
mkdir receptor
mkdir ligand
mkdir rstfiles

initialize_files () {
    cat > collect_$1.ptraj <<EOF
trajin ${trajin} 
$4 $5
trajout rstfiles/$1 restart multi
go
EOF


    cat > $2.sh <<EOF
#!/bin/bash 
#\$ -cwd
#
#\$ -m e
#\$ -S /bin/bash
#\$ -e error_file_ptraj_$1
#\$ -o output_file_ptraj_$1
#\$ -l h_rt=840:00:00
#\$ -q all.q
#\$ -pe orte 1

export LD_LIBRARY_PATH="/share/apps/pkgs/openmpi.1.8.1/lib:\${LD_LIBRARY_PATH}"
source /share/apps/pkgs/amber14/amber.sh
source /share/apps/pkgs/intel/bin/compilervars.sh intel64

echo Running on host `hostname`
echo Time is `date`
echo Directory is `pwd`

cpptraj ${prmtop_complex}.prmtop collect_$1.ptraj
EOF

    cat >$1.sh<<EOF
#!/bin/bash 
#\$ -cwd
#
#\$ -m e
#\$ -S /bin/bash
#\$ -e error_file_$1
#\$ -o output_file_$1
#\$ -l h_rt=840:00:00
#\$ -q all.q
#\$ -pe orte 1

export LD_LIBRARY_PATH="/share/apps/pkgs/openmpi.1.8.1/lib:\${LD_LIBRARY_PATH}"
source /share/apps/pkgs/amber14/amber.sh
source /share/apps/pkgs/intel/bin/compilervars.sh intel64

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

    qsub $2.sh
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
