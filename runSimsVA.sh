#!/bin/bash

#$ -q krt,krti,bio,pub64
#$ -pe openmp 64
#$ -ckpt restart
#$ -R y
##These jobs take A HECK OF A LOT of RAM, so avoid low-memory nodes:
#$ -l mem_size=512
module load krthornt/thorntonlab
cd $SGE_O_WORKDIR

model=$1
lambda=$2
outfile=$3
seed=$4
dominance=$5
sampler=$6

#We only do 32 sims at a time to (attempt to) keep peak RAM under control.  8*32 = 256 total replicates
/usr/bin/time -f "%e %M" -o $outfile.time python tennessen.py --model $model -l $lambda -o $outfile --seed $4 -d $5 --sampler $6 --batches 8 --cores 25 -t 100