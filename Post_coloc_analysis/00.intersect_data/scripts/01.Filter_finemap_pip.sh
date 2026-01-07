#!/bin/bash

for i in {1..22}
do

chr="chr${i}"

for file_num in $(seq -w 30000 10 30300); do
    
	formatted_file="$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results.txt"
    	
	if [ -f "$formatted_file" ]; then
		
		echo "Processing $formatted_file for chromosome ${i}"
	
		awk -F, '$1 > 0.01' ${formatted_file} > "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_pip.0.01.txt"
	
		echo "Creating bed file"

		cut -d$',' -f12 "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_pip.0.01.txt" | \
		awk -F$'\t' 'NR==FNR {vals[$1]; next} $4 in vals' - "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results.bed" \
		> "$home/stingseq_eqtl_overlap/results/UKBB_SuSiE_finemap/${file_num}/${file_num}_${chr}_finemap_results_pip.0.01.bed"

	else
		echo "File $formatted_file does not exist, skipping..."
	fi
done

done
