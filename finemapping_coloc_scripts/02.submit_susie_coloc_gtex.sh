#!/bin/bash
#SBATCH --job-name=susie_coloc_R                  # Job name
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -o ./Reports/coloc_out/coloc-%j.out # STDOUT
#SBATCH --array=1-22
#SBATCH --time=7-00:00:00


i=${SLURM_ARRAY_TASK_ID}

module purge
module load R/4.4.1

echo 'Running SuSiE COLOC Rscript'

NAME="GTEx"
EQTL="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/eQTL_sumstats/GTEx_whlbld.sumstats.chr${i}.hg19.txt"

# Create the directory if it doesn't already exist
mkdir -p "./Reports/coloc_out/${NAME}"

# Loop through the specified range of formatted.tsv files
for file_num in $(seq -w 30000 30300); do
    formatted_file="$home/stingseq_eqtl_overlap/data/UKBB_sumstats/${file_num}_formatted.tsv"
    if [ -f "$formatted_file" ]; then
        echo "Processing $formatted_file for chromosome ${i}"

        Rscript 02.susie_coloc_GTEx.R "${NAME}" "${i}" "$formatted_file" "${EQTL}"

    else
        echo "File $formatted_file does not exist, skipping..."
    fi
done

echo "Success for chromosome ${i}!"

# Move SLURM log files to the specified directory
mv "./Reports/coloc_out/coloc-${SLURM_JOB_ID}.out" "./Reports/coloc_out/${NAME}/"
