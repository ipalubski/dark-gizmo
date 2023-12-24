#!/bin/bash

## This script generates evrything needed to run a DM-only Isolated Halo 
## Positional Arguments:
## 1 Case Name
## 2 Cross Section
## 3 DM Force Softening
## 4 DM_g
## 5 DM_InteractionVelocityScale
## 6 Run Name (SLURM only)

##Create Directories
CASEDIR=$1
echo "$CASEDIR"
mkdir /lustre/home/ipalubski/simulations/$CASEDIR
mkdir /lustre/projects/palubski-group/${CASEDIR}out
OUTDIR=/lustre/projects/palubski-group/${CASEDIR}out
echo ${OUTDIR}

##Change Run Parameters
cp /lustre/home/ipalubski/dark-gizmo/base/* /lustre/home/ipalubski/simulations/$CASEDIR
cd /lustre/home/ipalubski/simulations/${CASEDIR}

sed -i "s/^DM_InteractionCrossSection\t.*/DM_InteractionCrossSection\t ${2}/" test.params
sed -i "s/^Softening_Type1\t.*/Softening_Type1\t ${3}/" test.params
sed -i "s/^DM_g\t.*/DM_g\t ${4}/" test.params
sed -i "s/^DM_InteractionVelocityScale\t.*/DM_InteractionVelocityScale\t ${5}/" test.params
sed -i "s@^%InitCondFile\t.*@%InitCondFile\t ${OUTDIR}/snapshot_001@" test.params
sed -i "s@^OutputDir\t.*@OutputDir\t ${OUTDIR}@" test.params

##Edit batch script
sed -i "s/SBATCH --job-name.*/SBATCH --job-name=${6}/" iso10.sub
sed -i "s@^OUTDIR.*@OUTDIR=${CASEDIR}out@" iso10.sub

##Compile Run
cd /lustre/home/ipalubski/dark-gizmo/
#make clean
#make
cp GIZMO /lustre/home/ipalubski/simulations/$CASEDIR/.

##Submit Job
#cd /lustre/home/ipalubski/simulations/$CASEDIR/
#sbatch iso10.sub
