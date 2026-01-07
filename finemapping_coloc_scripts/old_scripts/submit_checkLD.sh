#!/bin/bash
#SBATCH --job-name=CheckLDblocks                  # Job name
#--mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#--mail-user=sghatan@nygenome.org      # Where to send mail
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH -o ./Reports/CheckLD/checkld-%A-%a.out # STDOUT
#SBATCH --array=1-22
i=${SLURM_ARRAY_TASK_ID}

module purge
module load R/4.3.1

echo "Running  Rscript for chromosome ${i}"

echo "Running SuSiE FINEMAP Rscript for chromosome ${i}"

# Loop through the specified range of formatted.tsv files
for file_num in $(seq -w 30000 30001); do
    formatted_file="$home/stingseq_eqtl_overlap/data/UKBB_sumstats/${file_num}_formatted.tsv"
    if [ -f "$formatted_file" ]; then
        echo "Processing $formatted_file for chromosome ${i}"
        Rscript Check_LDblocks.R \
                "$formatted_file" \
                "$home/stingseq_eqtl_overlap/New_multi_coloc_pipeline/lead_loci/chr${i}.regions" \
                "$home/stingseq_eqtl_overlap/data/UKBB_LDmatrices/${i}/"
    else
        echo "File $formatted_file does not exist, skipping..."
    fi
done

echo "Success for chromosome ${i}!"
