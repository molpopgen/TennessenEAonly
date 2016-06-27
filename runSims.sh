#!/bin/bash

#$ -q krt,krti,bio,pub64
#$ -pe openmp 64

module load krthornt/thorntonlab
cd $SGE_O_WORKDIR

model=$1
lambda=$2
outfile=$3
seed=$4
dominance=$5
sampler=$6

#1024 replicates per parameter combo!!!
python tennessen.py --model $model -l $lambda -o $outfile --seed $4 -d $5 --sampler $6 --batches 16