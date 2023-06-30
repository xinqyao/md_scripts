#!/bin/bash
##SBATCH --account=chem36paid12
##SBATCH --account=CHEM36S12
#SBATCH --account=chem36r12
#SBATCH --partition=qGPU24
##SBATCH --partition=qCPU24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
##SBATCH --nodelist=acidsgcn001.rs.gsu.edu
##SBATCH --nodelist=acidsgcn001
#SBATCH --exclude=acidsgcn005,acidsgcn006
##SBATCH --exclude=acidsgn005
#SBATCH --gres=gpu:V100:1
##SBATCH --gres=gpu:RTX:1
##SBATCH --gres=gpu:Titan:1
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=xyao4@gsu.edu
#SBATCH --job-name=pin%ii%

echo
echo $SLURM_JOB_NODELIST
echo

bok=TRUE

nvidia-smi

module load CUDA/10.1.243-GCC-8.3.0 Amber20
#source /application/LocalBuilds/amber20_src/amber.sh 
#export LD_LIBRARY_PATH=/software/usr/cuda-8.0/lib64:/usr/lib64/openmpi/lib/:${LD_LIBRARY_PATH}
#export CUDA_VISIBLE_DEVICES="0"

###### Editable Parameters ######
rootdata=exchange
sys=ffpspr_pin1_box
nstlim=20000000

# ffpspr_pin1, ffpspr_ffpspr, ffpspr_tpp
atom1=2571
atom2=2582
atom3=2584
atom4=2593

# ffpspr_ppiase
#atom1=1819
#atom2=1830
#atom3=1832
#atom4=1841

# ffspr
#atom1=49
#atom2=57
#atom3=59
#atom4=68

# ffpspr
#atom1=49
#atom2=60
#atom3=62
#atom4=71
################################

#cd /scratch
mkdir $SLURM_JOB_ID
cd $SLURM_JOB_ID

if ! icd; then
   echo Error: irod command failed
   exit 1
fi
icd $rootdata
iget -r %dir%

## angles
angles=%angles%

# set rnum=1 but each run 40ns!!
for ((rnum=1; rnum<=1; rnum++)); do
#   for direc in 2.8 3.5 4.2 4.9 5.6; do
#   for direc in 3.5 4.2 4.9 5.6; do
   for direc in %dir%; do
#      if ! test -r $direc; then
#         mkdir $direc
#         tt=$( echo $direc | awk '{printf "%.8E", $1/4}' )
#         sed "s/8.33250000E+00/$tt/" ../prep/${sys}.prmtop > $direc/${sys}.prmtop
#      fi
      cd $direc

#      for ((angles=%a1%; angles>=%a2%; angles-=5)); do
      for ((iangle=%a1%; iangle<=%a2%; iangle++)); do
         f=${angles[iangle]}
         let rl=$f-180
         let rr=$f+180

         echo $rnum, $direc, $f $rl $rr

#         if ! test -r production.${f}.restrt; then
#            cp ../../windows/umbrella_${f}/production.A.restrt ./production.${f}.restrt
#         fi

         cat > noe.in <<EOF
# angle
 &rst
   ixpk= 0, nxpk= 0, iat=${atom1}, ${atom2}, ${atom3}, ${atom4}, r1=$rl, r2=$f, r3=$f, r4=$rr,
    rk2=32.83, rk3=32.83, ir6=1, ialtd=0,
 /
EOF

         if ! test -r production.in; then
         cat > production.in<<EOF
Initial minimization w/ position restraints on CypA, 9.0 cut
 &cntrl
  timlim = 999999, nmropt = 1, 
  ntx    = 5,       irest  = 1,       ntrx   = 1,      ntxo   = 2,
  ntpr   = 50000,      ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 1, ntwr = 100000,
  ntwprt = 0,

  ntf    = 2,       ntb    = 2,       igb  = 0,      dielc  = 1.0,
  cut    = 9.0,     nsnb   = 10,

  ipol   = 0,

  ibelly = 0,       ntr    = 0,

  imin   = 0,       
  maxcyc = 1000,   
  ncyc   = 5000,   
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = ${nstlim},
  nscm   = 10000, iwrap = 0,
  t      = 0.0,     dt     = 0.002,

  temp0  = 300.0,   tempi  = 100.0,   
  ig     = -1,   
  ntt    = 3,       
  tautp  = 2.0,  
  gamma_ln = 1.0,

  ntp    = 1,       pres0  = 1.0,     comp   = 44.6,
  taup   = 1.0,  
  barostat = 2,

  ntc    = 2,       tol    = 0.0001,

  jfastw = 0,       

  ivcap  = 0,       fcap   = 1.5,
 /
 &wt type='DUMPFREQ', istep1=10  /
 &wt type='END' /
DISANG=noe.in
DUMPAVE=ang.temp.out
EOF
         fi
 
         cp production.${f}.restrt production.run.restrt

         if ! pmemd.cuda -O -i production.in -p ${sys}.prmtop -c production.run.restrt -o production.${f}.out -r production.${f}.restrt; then
            bok=FALSE
         fi

         cat ang.temp.out >> ang.${f}.out
      done

#      if ! test -r wham; then 
#         mkdir -p wham/tmp
#         cp ../metafile wham/
#      fi
#      cd wham
#      for ifile in ../ang.[-0-9]*.out; do
#         awk '{a=$2; if(a<-60) a=a+360; if(a>300) a=a-360; print $1, a}' $ifile > tmp/$(basename $ifile)
#      done
#      /software/wham/wham/wham P -60 300 100 0.000001 300.0 0 metafile pmf.out 10 1234 &> err.log
#      awk 'NR > 1 && NR <102{print $1,$2}' pmf.out > fep.${rnum}.out
#      cd ..
      cd ..
   done
done

if test $bok == "FALSE"; then
   cd ..
   echo Done with errors.
   exit 1
fi

## check again incase the irod env changed
if ! icd; then
   echo Error: irod command failed
   echo Done with errors.
   exit 1
fi
icd $rootdata

if iput -rf %dir%; then
   cd ..
   rm -rf $SLURM_JOB_ID
   echo Done.
else
   cd ..
   echo Done with errors.
   exit 1
fi

