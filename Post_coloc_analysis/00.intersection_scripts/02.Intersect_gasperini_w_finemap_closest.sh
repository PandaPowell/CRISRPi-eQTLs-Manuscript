#!/bin/bash

# set dir
cd "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

module purge
module load bedtools
module load R/4.4.1

# Make temp folder for junk files
mkdir -p temp

# Make folder to hold credible sets
mkdir -p gasperini_credible_sets

# rm and initiate a file
echo -e "region\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site\tgwas" > gasperini_credible_sets/merged_gasperini_credset.txt

for i in {1..22}
do

chr="chr${i}"  # Assuming you have i set in your bash script.

awk -v chr="$chr" '{ if ($1 == chr) print $0 }' ../gasperini_spaceseq/resample_gRNA_target_sites.bed | \
sort -k1,1 -k2,2n > "temp/All_STING_seq_CRE_postions_${chr}.bed"

#2 We identify all the credible set of that variant by intersecting with finemapping results

for file_num in $(seq -w 30000 10 30300); do
    
	formatted_file="$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_pip.0.01.bed"
    	
	if [ -f "$formatted_file" ]; then
		
		mkdir -p gasperini_credible_sets/${file_num}		

		echo "Processing $formatted_file for chromosome ${i}"
	
		echo -e "region\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site" > "gasperini_credible_sets/${file_num}/${file_num}_${chr}_stingseq_credset.txt"

		awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' ${formatted_file} | sort -k1,1 -k2,2n > "temp/${file_num}_${chr}_finemap_results.bed"

		# Keep overlaps with distance less than 2000	
		bedtools closest -d -wb -a "temp/${file_num}_${chr}_finemap_results.bed" -b "temp/All_STING_seq_CRE_postions_${chr}.bed" | \
	 	awk '$10 < 2000' | sort -k12 | awk '!seen[$9]++' | sort -k5 > "temp/${file_num}_${chr}_overlap.bed"

		# Extract unique values from column 5 of File 2
		cut -d$'\t' -f5 "temp/${file_num}_${chr}_overlap.bed" | sort -u > temp/${file_num}_${chr}_credible_set.txt

		# Filter File 1 to include only lines where the value in column 5 matches those in unique_values.txt
		awk -F$'\t' 'NR==FNR {vals[$1]; next} $5 in vals' "temp/${file_num}_${chr}_credible_set.txt" "temp/${file_num}_${chr}_finemap_results.bed" | \
		sort -k5 > "temp/${file_num}_${chr}_finemap_filtered.txt"

		join -t $'\t' -1 5 -2 5 -o auto "temp/${file_num}_${chr}_finemap_filtered.txt" "temp/${file_num}_${chr}_overlap.bed" | \
                awk -v OFS='\t' '{print $1,$2,$5,$9,$11,$12,$13}' >> "gasperini_credible_sets/${file_num}/${file_num}_${chr}_stingseq_credset.txt"

		join -t $'\t' -1 5 -2 5 -o auto "temp/${file_num}_${chr}_finemap_filtered.txt" "temp/${file_num}_${chr}_overlap.bed" | \
		awk -v OFS='\t' -v gwas="$file_num" '{print $1,$2,$5,$9,$11,$12,$13,gwas}' >> "gasperini_credible_sets/merged_gasperini_credset.txt"

		
	fi
done

done

rm -r temp

# header of final output file will be: 
# region,chr,finemap_snp_lower,finemap_snp_upper,finemap_snp, - finemapped snps from full credible set
# chr,finemap_snp_lower,finemap_snp_upper,finemap_snp, - fine mapped snp that originally intersected grna region
# chr,lower_grna_target,upper_grna_target,region_gene_pair - gasperini regions
