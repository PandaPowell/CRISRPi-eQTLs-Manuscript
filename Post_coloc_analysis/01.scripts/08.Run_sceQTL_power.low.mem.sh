#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --array=0-3  # Adjust the range to match the number of elements in your array
#SBATCH --output=slurm-%A_%a.out  # Save output for each array job
#SBATCH --mem=50G

module purge
module load R/4.4.1

# Define the array
names=("B" "Mono" "NK" "DC")

# Get the current array task's index
name=${names[$SLURM_ARRAY_TASK_ID]}

# Run the R script
Rscript 08.calculate_sceqtl_power.R ${name}

echo "Done!"
