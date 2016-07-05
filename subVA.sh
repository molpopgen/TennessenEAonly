#!/bin/bash

sampler=VA
#Remember: GBR model = dominance doesn't matter



seed=$RANDOM

for model in additive #gbr
do
    for lambda in 0.025 0.05 0.1 0.25 0.5
    do
	ofile=$model."lambda"$lambda.$sampler.h5
	name=$model$sampler
	qsub -N $name runSimsVA.sh $model $lambda $ofile $seed 1.0 $sampler
	seed=$RANDOM
    done
done

seed=$RANDOM
for model in multi
do
    for lambda in 0.025 0.05 0.1 0.25 0.5
    do
	for dominance in 0 0.25 1
	do
	    ofile=$model."lambda"$lambda.$sampler.h$dominance.h5
	    name=$model$sampler
	    qsub -N $name runSimsVA.sh $model $lambda $ofile $seed $dominance $sampler
	    seed=$RANDOM
	done
    done
done