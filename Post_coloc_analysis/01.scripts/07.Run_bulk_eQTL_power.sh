#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00  # 7 days

module purge
module load R/4.4.1

# Run the R script
Rscript 08.calculate_bulk_eqtl_power.R

echo "Done!"
