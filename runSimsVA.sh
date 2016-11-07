#!/bin/bash

#$ -q krt,krti,bio,pub64,abio,free64
#$ -pe openmp 4
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
/usr/bin/time -f "%e %M" -o $outfile.$SGE_TASK_ID.time python -u tennessen.py --model $model -l $lambda -o $outfile.$SGE_TASK_ID.h5 --seed $SEED -d $5 --sampler $6 --batches 1 --cores 1 -t 1000 --bigstub $outfile.$SGE_TASK_ID
