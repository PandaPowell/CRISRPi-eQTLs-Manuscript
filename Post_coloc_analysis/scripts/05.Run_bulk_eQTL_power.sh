#!/bin/bash
#SBATCH --job-name=eqtl_power
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00           # 7 days
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail
mkdir -p logs

module purge
module load R/4.4.1

# Run the R scripts (sequentially)
Rscript power_scripts/08.calculate_bulk_eqtl_power.R
Rscript power_scripts/08.calculate_interval_power.R

echo "Done!"
