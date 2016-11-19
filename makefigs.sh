#!/bin/bash

#$ -q krt,krt,bio 
#$ -pe openmp 4

cd $SGE_O_WORKDIR
module load krthornt/thorntonlab
ls *.h5|grep summ|xargs rm -f

jupyter nbconvert --execute --to notebook --ExecutePreprocessor.timeout=6000 Stats.ipynb
