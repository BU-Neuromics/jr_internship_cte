#!/bin/bash
#$ -P cteseq
#$ -pe omp 32
#$ -j y

cd /restricted/projectnb/cteseq/jrose
module load R
Rscript CTE_deseq_cont_group_qsub_high.R

