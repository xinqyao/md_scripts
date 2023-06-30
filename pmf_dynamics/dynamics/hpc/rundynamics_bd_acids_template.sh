#!/bin/bash
##SBATCH --account=chem36paid12
##SBATCH --account=CHEM36S12
#SBATCH --account=chem36r12
##SBATCH --partition=qGPU24
#SBATCH --partition=qCPU24
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
####SBATCH --nodelist=acidsgcn001.rs.gsu.edu
##SBATCH --nodelist=acidsgcn001
##SBATCH --exclude=acidsgcn005
##SBATCH --gres=gpu:V100:1
#SBATCH --gres=gpu:RTX:1
##SBATCH --gres=gpu:Titan:1
##SBATCH --mail-type=END,FAIL
##SBATCH --mail-user=xyao4@gsu.edu
#SBATCH --job-name=pin%ii%

#export LD_LIBRARY_PATH=/software/usr/cuda-8.0/lib64:/usr/lib64/openmpi/lib/:${LD_LIBRARY_PATH}
#export CUDA_VISIBLE_DEVICES="0"

echo
echo $SLURM_JOB_NODELIST
echo

bok=TRUE

nvidia-smi

module load CUDA/10.1.243-GCC-8.3.0 Amber20

#trap "exit" INT TERM
#trap "kill 0" EXIT

#cd /scratch
mkdir $SLURM_JOB_ID
cd $SLURM_JOB_ID

if ! icd; then
   echo irod command failed
   exit 1
fi
icd exchange
iget -r %dir%

prv=180
#set simult = (2500000 5000000 10000000 25000000 37500000 50000000)
#simult=(1000000 2500000 5000000 10000000 50000000)
#tnum=4

prmtop_in=%prmtop_in%
prmtop=%prmtop%

#for direc in 2.8 3.5 4.2 4.9 5.6; do
for direc in %dir%; do
#  if ! test -r $direc; then
#     mkdir $direc
#  fi

  cd $direc

#  echo $direc, $tnum, ${simult[$tnum]}
  echo $direc

#  if ! test -r production.${prv}.restrt; then
#     cp ../../windows/umbrella_${prv}/production.A.restrt ./production.${prv}.restrt
#  fi

#  if ! test -r ffpspr.bd.prmtop; then 
#    tt=$( echo 14.0 | awk '{printf "%.8E", $1/4}' )
#    sed "s/8.33250000E+00/$tt/" ../../prep/ffpspr_ppiase_box.prmtop > ./ffpspr.bd.prmtop
#  fi

#  if ! test -r ffpspr_box.prmtop; then 
#    tt=$( echo $direc | awk '{printf "%.8E", $1/4}' )
#    sed "s/8.33250000E+00/$tt/" ../../prep/ffpspr_ppiase_box.prmtop > ./ffpspr_box.prmtop
#  fi
  
  cat > production.bd.in<<eofxx
Initial minimization w/ position restraints on CypA, 9.0 cut
 &cntrl
  timlim = 999999, nmropt = 0, 
  ntx    = 5,       irest  = 0,       ntrx   = 1,      ntxo   = 2,
  ntpr   = 50000,      ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
  ioutfm = 1, ntwr = -500000,
  ntwprt = 0,

  ntf    = 2,       ntb    = 2,       igb  = 0,      dielc  = 1.0,
  cut    = 9.0,     nsnb   = 10,

  ipol   = 0,

  ibelly = 0,       ntr    = 0,

  imin   = 0,       
  maxcyc = 1000,   
  ncyc   = 5000,   
  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,

  nstlim = 50000000,
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
eofxx

  cat > noe.in <<eof
# angle
 &rst
   ixpk= 0, nxpk= 0, iat=1819, 1830, 1832, 1841, r1=0, r2=180, r3=180, r4=360,
    rk2=0.0, rk3=0.0, ir6=1, ialtd=0,
 /
eof

# set number=1 but run 100ns; save restrt file every 1ns!
  for ((number=1; number<=1; number++)); do

#     echo $direc, $tnum, ${simult[$tnum]}, $number
     echo $direc, $number
     
     if test ! -d restrt; then
        mkdir restrt
     fi

     if ! pmemd.cuda -O -i production.bd.in -p ${prmtop}.bd.prmtop -c production.${prv}.restrt -o production.bd.out -r restrt/production.bd.restrt; then
        bok=FALSE
     fi

     ## convert restrt file names ##
     cd restrt
     j=100
     for i in `ls -t production.bd.restrt_*`; do 
        #echo $i
        mv $i production.bd.${j}.restrt
        let j=j-1
     done
     cd ..

#     sleep 15

     cp restrt/production.bd.restrt production.${prv}.restrt

#     cat > production.in<<eofx
#Initial minimization w/ position restraints on CypA, 9.0 cut
# &cntrl
#  timlim = 999999, nmropt = 1, 
#  ntx    = 5,       irest  = 1,       ntrx   = 1,      ntxo   = 1,
#  ntpr   = 50000,      ntwx   = 0,       ntwv   = 0,      ntwe   = 0,
#  ioutfm = 0, ntwr = 50000,
#  ntwprt = 0,
#
#  ntf    = 2,       ntb    = 2,       igb  = 0,      dielc  = 1.0,
#  cut    = 9.0,     nsnb   = 10,
#
#  ipol   = 0,
#
#  ibelly = 0,       ntr    = 0,
#
#  imin   = 0,       
#  maxcyc = 1000,   
#  ncyc   = 5000,   
#  ntmin  = 1,       dx0    = 0.1,     drms   = 0.0001,
#
#  nstlim = ${simult[tnum]},
#  nscm   = 10000, iwrap = 1,
#  t      = 0.0,     dt     = 0.002,
#
#  temp0  = 300.0,   tempi  = 100.0,   
#  ig     = -1,   
#  ntt    = 3,       
#  tautp  = 2.0,  
#  gamma_ln = 1.0,
#
#  ntp    = 1,       pres0  = 1.0,     comp   = 44.6,
#  taup   = 1.0,  
#  barostat = 2,
#
#  ntc    = 2,       tol    = 0.0001,
#
#  jfastw = 0,       
#
#  ivcap  = 0,       fcap   = 1.5,
# /
# &wt type='DUMPFREQ', istep1=10  /
# &wt type='END'  /
#DISANG=noe.in
#DUMPAVE=ang.${number}.out
#eofx

#     pmemd.cuda -O -i production.in -p ffpspr_box.prmtop -c production.bd.restrt -o production.out -r production.restrt &
#     pmemd.cuda -O -i production.in -p ffpspr_box.prmtop -c production.bd.restrt -o production.out -r production.restrt 

#     sleep 60

#     while test -r ../isodetect; do
#        if test -r iso.${number}.out.pre; then
#           if ../isodetect ang.${number}.out iso.${number}.out 175 100 9 iso.${number}.out.pre; then
#              break
#           fi
#        else
#           if ../isodetect ang.${number}.out iso.${number}.out 175 100 9; then
#              break
#           fi
#        fi
#        mv iso.${number}.out iso.${number}.out.pre
#        sleep 60
#     done

#     if test $( jobs -pr | wc -l ) -gt 0; then
#        kill -15 $(jobs -pr)
#        wait
#     fi

#     sleep 15

  done

  echo "Exiting directory ..."

  cd ..

#@ tnum = $tnum + 1
#let tnum=tnum+1

done

if test $bok == "FALSE"; then
   cd ..
   echo Done with errors.
   exit 1
fi

if iput -rf %dir%; then
   cd ..
   rm -rf $SLURM_JOB_ID
   echo Done.
else
   cd ..
   echo Done with errors.
   exit 1
fi
