#!/bin/bash
#SBATCH --job-name=susie_finemap_R                  # Job name
#SBATCH --mem=50G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH -o ./Reports/finemap_out/finemap-%A-%a.out # STDOUT
#SBATCH --array=13

i=${SLURM_ARRAY_TASK_ID}

# chr1-22 50G

REGION_FILE="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/gwas_regions/merged_blood_trait_regions.txt"

module purge
module load R/4.4.1

echo "Running SuSiE FINEMAP Rscript for chromosome ${i}"

# Loop through the specified range of formatted.tsv files
for file_num in $(seq -w 30000 30300); do
    formatted_file="$home/stingseq_eqtl_overlap/data/UKBB_sumstats/${file_num}_formatted.tsv"
    if [ -f "$formatted_file" ]; then
        echo "Processing $formatted_file for chromosome ${i}"

	mkdir -p "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}"

        Rscript 01.susie_finemap.R \
                "$formatted_file" \
                "$home/stingseq_eqtl_overlap/results/gwas_regions/merged_blood_trait_regions.txt" \
                "$home/stingseq_eqtl_overlap/data/UKBB_LDmatrices/" \
		${i}
    else
        echo "File $formatted_file does not exist, skipping..."
    fi
done

echo "Success for chromosome ${i}!"
