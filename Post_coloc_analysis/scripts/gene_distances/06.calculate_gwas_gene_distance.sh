#!/bin/bash

module load bedtools

sort -k1,1 -k2,2n "../data/gencode.v33lift37.annotation.bed" > "../data/gencode.v33lift37.annotation.sorted.bed"

for i in {1..22}
do

chr="chr${i}"

	for file_num in $(seq -w 30000 10 30300); do

        	formatted_file="$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_pip.0.01.txt"

        	if [ -f "$formatted_file" ]; then
		
		awk -F, '{print"chr"$13"\t"$17-1"\t"$17"\t"$12"\t"$37"\t"$1}' ${formatted_file} | tail -n +2 | sort -k6,6nr | awk '!seen[$5]++' | sort -k1,1 -k2,2n \
		> "../../results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_leadpip.filtered.sorted.bed" 

		bedtools closest -D a -wb -a "../../results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_leadpip.filtered.sorted.bed" \
		-b "../data/gencode.v33lift37.annotation.sorted.bed" \
		> "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_gene_distance.bed"
		
		fi
	done

done

# Define output path
merged_dir="$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/merged"
mkdir -p "$merged_dir"

# Define header line
header="chr\tpos_l\tpos\tsnp\tregion_cs\tpip\tgene_chr\ttss_l\ttss\tgene\tdistance"

# Write header and concatenate files
{
  echo -e "$header"
  find "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap" \
    -name "*_gene_distance.bed" \
    -exec cat {} +
} > "$merged_dir/all_gene_distance.bed"
