#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module purge
module load R/4.3.1

#Rscript download_data.R

Rscript download_sumstats.R
