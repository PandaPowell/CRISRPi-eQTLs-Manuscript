#!/bin/bash

# set dir
cd "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

module load bedtools

# Make temp folder for junk files
mkdir -p temp

# Make folder to hold credible sets
mkdir -p sting_seq_credible_sets/BCX
mkdir -p gasperini_credible_sets/BCX
mkdir -p encode_credible_sets/BCX

array=("BAS" "EOS" "HCT" "HGB" "LYM" "MCHC" "MCH" "MCV" "MON" "MPV" "NEU" "PLT" "RBC" "RDW" "WBC")

sort -k1,1 -k2,2n ../CRISPR_data/All_Stingseq_pos.bed > "temp/All_STING_seq_CRE_postions.bed"

sort -k1,1 -k2,2n ../gasperini_spaceseq/resample_gRNA_target_sites.bed > "temp/All_gasperini_CRE_postions.bed"

awk '{print $2"\t"$3"\t"$4"\t"$1}' ../CRISPR_data/NoGasperini_crispri_data.tsv | tail -n +2 | \
sort -k1,1 -k2,2n > temp/All_encode_CRE_positions.bed

file_array=("temp/All_STING_seq_CRE_postions.bed" "temp/All_gasperini_CRE_postions.bed"  "temp/All_encode_CRE_positions.bed")

output_array=("sting_seq_credible_sets/BCX/" "gasperini_credible_sets/BCX/" "encode_credible_sets/BCX/")

echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tSNP_coord\tgwas" > "sting_seq_credible_sets/BCX/merged_credset.txt"

echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site\tgwas" \
> "gasperini_credible_sets/BCX/merged_credset.txt"

echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site\tgwas" \
> "encode_credible_sets/BCX/merged_credset.txt"


for i in "${array[@]}"; do
    
	formatted_file="$home/stingseq_eqtl_overlap/results/MAGE_finemap/${i}_finemap_results.bed"
    	
	if [ -f "$formatted_file" ]; then
		
		echo "Processing $formatted_file for ${i}"
	
		awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $4, $5}' ${formatted_file} | sort -k1,1 -k2,2n > "temp/${i}_finemap_results.bed"
	
		echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tSNP_coord\tgwas" > "sting_seq_credible_sets/BCX/${i}_credset.txt"

                echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site\tgwas" > "gasperini_credible_sets/BCX/${i}_credset.txt"

		echo -e "Sentinel_MarkerName\tchr\tfinemap_snp\tfinemap_snp_intersect_grna\tlower_grna_target_site\tupper_grna_target_site\ttarget_site\tgwas" > "encode_credible_sets/BCX/${i}_credset.txt"

		for j in {0..2}; do

			echo ${file_array[$j]}

			bedtools intersect -wb -a "temp/${i}_finemap_results.bed" -b ${file_array[$j]} | \
	 		sort -k5 > "temp/${i}_overlap.bed"

			# Extract unique values from column 5 of File 2
			cut -d$'\t' -f5 "temp/${i}_overlap.bed" | sort -u > temp/${i}_credible_set.txt

			# Filter File 1 to include only lines where the value in column 5 matches those in unique_values.txt
			awk -F$'\t' 'NR==FNR {vals[$1]; next} $5 in vals' "temp/${i}_credible_set.txt" "temp/${i}_finemap_results.bed" | \
			sort -k5 > "temp/${i}_finemap_filtered.txt"

			if [ ${file_array[$j]} = "temp/All_STING_seq_CRE_postions.bed" ]; then

				join -t $'\t' -1 5 -2 5 -o auto "temp/${i}_finemap_filtered.txt" "temp/${i}_overlap.bed" | \
				awk -v OFS='\t' -v gwas="$i" '{print $1,$2,$5,$9,$13,gwas}' \
				>> "${output_array[$j]}${i}_credset.txt"

				join -t $'\t' -1 5 -2 5 -o auto "temp/${i}_finemap_filtered.txt" "temp/${i}_overlap.bed" | \
                                awk -v OFS='\t' -v gwas="$i" '{print $1,$2,$5,$9,$13,gwas}' \
                                >> "${output_array[$j]}merged_credset.txt"

			else

				join -t $'\t' -1 5 -2 5 -o auto "temp/${i}_finemap_filtered.txt" "temp/${i}_overlap.bed" | \
	                	awk -v OFS='\t' -v gwas="$i" '{print $1,$2,$5,$9,$11,$12,$13,gwas}' \
				>> "${output_array[$j]}${i}_credset.txt"

				join -t $'\t' -1 5 -2 5 -o auto "temp/${i}_finemap_filtered.txt" "temp/${i}_overlap.bed" | \
                                awk -v OFS='\t' -v gwas="$i" '{print $1,$2,$5,$9,$11,$12,$13,gwas}' \
                                >> "${output_array[$j]}merged_credset.txt"

			fi
		done
	fi
done

rm -r temp
