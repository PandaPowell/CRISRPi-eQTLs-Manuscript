#!/bin/bash
#SBATCH --mem=40G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --array=1-22

i=${SLURM_ARRAY_TASK_ID}

module purge
module load bedtools

mkdir -p bedfiles
mkdir -p UKBB_susie_finemap_hg38

DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap_lbf/"

for file_num in $(seq -w 30000 10 30300); do

  formatted_file="${DIR}${file_num}/${file_num}_chr${i}_finemap_results.txt"

  if [ -f "$formatted_file" ]; then

  mkdir -p UKBB_susie_finemap_hg38/${file_num}

  bed_file="${DIR}${file_num}/${file_num}_chr${i}_finemap_results.bed"

  echo "${bed_file}"

  # First step is to obtain rsid from position ranges
  bedtools intersect -wa -wb -a "${bed_file}" -b $home/genome/dbsnp151_GRCh37p13/sorted_bed_chr_${i}.bed.gz > bedfiles/hg19_rsids_chr${i}.txt
 
  # Use vlookup to find rsids and pos in hg38 files
  awk -v FS="\t" 'FNR==NR {a[$8] = $0; next} $4 in a {print a[$4], $0}' bedfiles/hg19_rsids_chr${i}.txt <(zcat $home/genome/dbsnp151_GRCh38p7/bed_chr_${i}.bed.gz) | 
  awk 'BEGIN {OFS=","} {print$1,$3,$4,$8,$13}' > bedfiles/hg19_hg38_rsids_chr${i}.txt

  echo -e "chr,pos_hg19,variant,rsid,hg38_pos,variant.x,lbf_1,lbf_2,lbf_3,lbf_4,lbf_5,lbf_6,lbf_7,lbf_8,lbf_9,lbf_10,chr.x,region,chr.y,pos,a0,a1,pval" \
  > "UKBB_susie_finemap_hg38/${file_num}/${file_num}.chr${i}_finemap_results.hg38.txt"

  # Use vlookup to join rsids file back to sumstats file based on SNP column and filter columns of finemap file
  awk 'BEGIN { OFS="," } 
  FNR==NR {FS=","; a[$3] = $0; next} 
  FNR!=NR {FS=","; if ($1 in a) {print a[$1], $0}}' \
  bedfiles/hg19_hg38_rsids_chr${i}.txt ${formatted_file} \
  >> "UKBB_susie_finemap_hg38/${file_num}/${file_num}.chr${i}_finemap_results.hg38.txt"

  fi

done

echo "Done!"
