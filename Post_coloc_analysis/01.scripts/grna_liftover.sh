#!/bin/bash
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

module load bedtools

mkdir -p temp

for i in {1..22}; do

  echo "chr${i}"

  # Fixing the awk variable expansion with -v
  tail -n +2 cres_with_grnas.txt | awk -F, -v i="$i" 'BEGIN {OFS="\t"} $11 == i {print "chr"$11,$10-1,$10,$12}' | sort | uniq > temp/temp_chr${i}.bed

  # First step is to obtain rsid from position ranges
  bedtools intersect -wa -wb -a temp/temp_chr${i}.bed -b $home/genome/dbsnp151_GRCh37p13/sorted_bed_chr_${i}.bed.gz > temp/hg19_rsids_chr${i}.txt

  # Use vlookup to find rsids and pos in hg38 files
  awk -v FS="\t" 'FNR==NR {a[$8] = $0; next} $4 in a {print a[$4], $0}' temp/hg19_rsids_chr${i}.txt <(zcat $home/genome/dbsnp151_GRCh38p7/bed_chr_${i}.bed.gz) |
  awk 'BEGIN {OFS=","} {print $1,$3,$4,$8,$13}' > temp/hg19_hg38_rsids_chr${i}.txt

done
