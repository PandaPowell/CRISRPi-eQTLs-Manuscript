#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --array=30000-30300:10

# Then, in your script:
i=${SLURM_ARRAY_TASK_ID}

module purge
module load R/4.4.1

echo 'Running SuSiE COLOC Rscript'

Rscript 01.eQTL_catalogue.R ${i}
