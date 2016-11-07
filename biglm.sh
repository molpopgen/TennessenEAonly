#!/bin/bash

#$ -q krt,krti,pub64,bio,free64
#$ -pe openmp 3
#$ -t 1-1500
#$ -ckpt restart

module load R

cd $SGE_O_WORKDIR

INFILE=`head -n $SGE_TASK_ID bigN/infiles | tail -n 1`
OUTFILEBASE=`basename $INFILE .gz`

Rscript va.R $INFILE bigN/$OUTFILEBASE.VA.gz


