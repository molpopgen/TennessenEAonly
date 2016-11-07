#!/bin/bash

#$ -q krt,krti,bio,pub64,abio,free64
#$ -pe openmp 1
#$ -ckpt restart
#$ -R y
##These jobs take A HECK OF A LOT of RAM, so avoid low-memory nodes:
#$ -l mem_size=512
#$ -t 1-100

module load krthornt/thorntonlab
cd $SGE_O_WORKDIR

model=$1
lambda=$2
outfile=$3
seed=$4
dominance=$5
sampler=$6

#1 rep at a time, but use an array job
SEED=`echo "$RANDOM*$SGE_TASK_ID"|bc -l`
python tennessen.py --model $model -l $lambda -o $outfile.$SGE_TASK_ID.h5 --seed $SEED -d $5 --sampler $6 --batches 1 --cores 1 -t 1000 --bigstub bigN/$outfile.$SGE_TASK_ID
