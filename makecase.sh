#!/bin/bash

## This script generates evrything needed to run a DM-only Isolated Halo 
## Positional Arguments:
## 1 Case Name
## 2 Cross Section
## 3 DM Force Softening
## 4 GRAVHFactor
## 5 SIDM SmoothingFactor
## 6 Run Name (SLURM only)

##Create Directories
CASEDIR=$1
echo "$CASEDIR"
mkdir /pub/ipalubsk/gizmo-arepo/$CASEDIR
mkdir /pub/ipalubsk/${CASEDIR}out
OUTDIR=/pub/ipalubsk/${CASEDIR}out
echo ${OUTDIR}

##Change Run Parameters
cp /pub/ipalubsk/gizmo-arepo/base/* /pub/ipalubsk/gizmo-arepo/$CASEDIR
cd /pub/ipalubsk/gizmo-arepo/${CASEDIR}

sed -i "s/^DM_InteractionCrossSection\t.*/DM_InteractionCrossSection\t ${2}/" test.params
sed -i "s/^Softening_Type1\t.*/Softening_Type1\t ${3}/" test.params
sed -i "s/^GRAVHfactor\t.*/GRAVHfactor\t ${4}/" test.params
sed -i "s/^SIDMSmoothingFactor\t.*/SIDMSmoothingFactor\t ${5}/" test.params
sed -i "s@^OutputDir\t.*@OutputDir\t ${OUTDIR}@" test.params

##Edit batch script
sed -i "s/SBATCH --job-name.*/SBATCH --job-name=${6}/" iso10.sub
sed -i "s@^OUTDIR.*@OUTDIR=${CASEDIR}out@" iso10.sub

##Compile Run
cd /pub/ipalubsk/gizmo-arepo/
#make clean
#make
cp GIZMO /pub/ipalubsk/gizmo-arepo/$CASEDIR/.

##Submit Job
#cd /pub/ipalubsk/gizmo-public/$CASEDIR/
#sbatch iso10.sub
