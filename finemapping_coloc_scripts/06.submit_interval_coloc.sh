#!/bin/bash
#SBATCH --job-name=susie_coloc_R                  # Job name
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -o ./Reports/coloc_out/Interval_coloc-%j.out # STDOUT
#SBATCH --array=1-22
#SBATCH --time=7-00:00:00


i=${SLURM_ARRAY_TASK_ID}

module purge
module load R/4.4.1

echo 'Running SuSiE COLOC Rscript'

NAME="Interval"
EQTL="/gpfs/commons/datasets/controlled/INTERVAL/public_sumstats/cis/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_nominal_chr${i}.tsv"

# Create the directory if it doesn't already exist
mkdir -p "./Reports/coloc_out/${NAME}"

# Loop through the specified range of formatted.tsv files
for file_num in $(seq -w 30000 30300); do
    formatted_file="$home/stingseq_eqtl_overlap/data/UKBB_sumstats/${file_num}_formatted.tsv"
    if [ -f "$formatted_file" ]; then
        echo "Processing $formatted_file for chromosome ${i}"

        Rscript 06.Interval_coloc.R "${NAME}" "${i}" "$formatted_file" "${EQTL}"

    else
        echo "File $formatted_file does not exist, skipping..."
    fi
done

echo "Success for chromosome ${i}!"

# Move SLURM log files to the specified directory
mv "./Reports/coloc_out/coloc-${SLURM_JOB_ID}.out" "./Reports/coloc_out/${NAME}/"
