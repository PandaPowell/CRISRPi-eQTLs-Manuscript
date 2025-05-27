#!/bin/bash
#SBATCH --job-name=parquet_filter
#SBATCH --output=parquet_filter_%j.out  # Output log
#SBATCH --error=parquet_filter_%j.err   # Error log
#SBATCH --cpus-per-task=16               # Adjust CPU cores as needed
#SBATCH --mem=40G                        # Adjust memory as needed
#SBATCH --time=7-00:00:00  # 7 days

module load python

# Convert Jupyter notebook to Python script using micromamba
/gpfs/commons/home/sghatan/.local/bin/micromamba run -n crisprQTL jupyter nbconvert --to script Extract_MetaLCL.ipynb

# Run the converted Python script using micromamba
/gpfs/commons/home/sghatan/.local/bin/micromamba run -n crisprQTL python Extract_MetaLCL.py

echo "Done"
