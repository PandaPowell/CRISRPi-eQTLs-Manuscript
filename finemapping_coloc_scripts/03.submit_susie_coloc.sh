#!/bin/bash
#SBATCH --job-name=susie_coloc_R                  # Job name
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -o ./Reports/coloc_out/coloc-%j.out # STDOUT
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}

module purge
module load R/4.4.1

echo 'Running SuSiE COLOC Rscript'

cell_types=$1
path=$2
BASE="/gpfs/commons/groups/lappalainen_lab/woliveros/231005_OneK1K/analysis/02_QTL_calling/lvl1Res/"

NAME="${cell_types}"
EQTL="${BASE}${path}.chr${i}.parquet"

# Create the directory if it doesn't already exist
mkdir -p "./Reports/coloc_out/${NAME}"

echo "${cell_types} : ${EQTL}"

# Loop through the specified range of formatted.tsv files
for file_num in $(seq -w 30000 30300); do
    formatted_file="$home/stingseq_eqtl_overlap/data/UKBB_sumstats/${file_num}_formatted.tsv"
    if [ -f "$formatted_file" ]; then
        echo "Processing $formatted_file for chromosome ${i}"

        Rscript susie_coloc_beta.2.1.R "${NAME}" "${i}" "$formatted_file" "${EQTL}"

    else
        echo "File $formatted_file does not exist, skipping..."
    fi
done

echo "Success for chromosome ${i}!"

# Move SLURM log files to the specified directory
mv "./Reports/coloc_out/coloc-${SLURM_JOB_ID}.out" "./Reports/coloc_out/${NAME}/"
