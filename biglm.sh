#!/bin/bash

#$ -q krt,krti,pub64,sf,bio
#$ -pe openmp 3
#$ -t 1-1500

module load R

cd $SGE_O_WORKDIR

INFILE=`head -n $SGE_TASK_ID bigN/infiles | tail -n 1`
OUTFILEBASE=`basename $INFILE .gz`

Rscript va.R $INFILE bigN/$OUTFILEBASE.VA.gz


