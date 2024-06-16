#!/bin/bash
#$ -P cteseq
#$ -pe omp 32
#$ -j y
#$ -V
 
cd /restricted/projectnb/cteseq/jrose
module load R
Rscript CTE_deseq_ARTB_2021_model_qsub.R

