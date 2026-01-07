#!/bin/bash

# set dir
cd "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

mkdir -p temp

echo "GWAS_variant,region,chr,finemap_snp_intersect_grna,SS_coord,gwas,nsnps,gwas_hit_hg38,eqtl_hit_hg38,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,idx1,idx2,eqtl_name,\
gwas_name,study_id,dataset_id,molecular_id,gwas_pval,eqtl_hit_hg19" \
> "stingseq_eqtl_overlap/eqtl_catalogue_coloc_stingseq.txt"

echo "GWAS_variant,region,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,gwas_num,nsnps,gwas_hit_hg38,eqtl_hit_hg38,PP.H0.abf,\
PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,idx1,idx2,eqtl_name,gwas_name,study_id,dataset_id,molecular_id,gwas_pval,eqtl_hit_hg19" \
> "gasperini_eqtl_overlap/eqtl_catalogue_coloc_gasperini.txt"

echo "GWAS_variant,region,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,gwas_num,nsnps,gwas_hit_hg38,eqtl_hit_hg38,PP.H0.abf,\
PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,idx1,idx2,eqtl_name,gwas_name,study_id,dataset_id,molecular_id,gwas_pval,eqtl_hit_hg19" \
> "encode_eqtl_overlap/eqtl_catalogue_coloc_crispri.txt"

for file_num in $(seq -w 30000 10 30300); do

	file_array=("sting_seq_credible_sets/merged_stingseq_credset.txt" "gasperini_credible_sets/merged_gasperini_credset.txt" "encode_credible_sets/merged_encode_credset.txt")
	output_array=("stingseq_eqtl_overlap/eqtl_catalogue_coloc_stingseq.txt" "gasperini_eqtl_overlap/eqtl_catalogue_coloc_gasperini.txt" \
"encode_eqtl_overlap/eqtl_catalogue_coloc_crispri.txt")

	echo "Running for $file_num"

	for i in "${!file_array[@]}"; do

		coloc_file="$home/stingseq_eqtl_overlap/eQTL_catalogue/coloc_results/${file_num}_coloc_results.txt"
		cred_set_file="${file_array[$i]}"
		
		echo "Running for ${file_array[$i]}"

                # Extract unique SNP ids from column 1 of coloc file
                cut -d',' -f16 $coloc_file | \
                sort -u > temp/${file_num}_coloc_snps.txt

		# Sort coloc file
                sort -t',' -k16 $coloc_file > "temp/catalogue_${file_num}_coloc_results.txt"
		
		if [ "$cred_set_file" == "sting_seq_credible_sets/merged_stingseq_credset.txt" ]; then
			# Filter stingseq credible set file to include only lines where SNP from coloc file is present and gwas number
			awk -v OFS=',' -v file_num="$file_num" 'NR==FNR {vals[$1]; next} $3 in vals && $6 == file_num \
			{print $1, $2, $3, $4, $5, $6}' "temp/${file_num}_coloc_snps.txt" "${cred_set_file}" | \
			sort -t',' -k3 > "temp/${file_num}_finemap_filtered.txt"

                	# Check if there are any matching SNP ids if so then we join the files
                	if grep -qf temp/${file_num}_coloc_snps.txt <(cut -d',' -f3 "temp/${file_num}_finemap_filtered.txt"); then

                        	# Perform join if matches are found

                        	join -t ',' -1 3 -2 16 -o auto "temp/${file_num}_finemap_filtered.txt" \
                        	"temp/catalogue_${file_num}_coloc_results.txt" >> "${output_array[$i]}"
			fi

		else 

			# Filter stingseq credible set file to include only lines where SNP from coloc file is present
                        awk -v OFS=',' -v file_num="$file_num" 'NR==FNR {vals[$1]; next} $3 in vals && $8 == file_num \
			{print $1, $2, $3, $4, $5, $6, $7, $8}' "temp/${file_num}_coloc_snps.txt" $cred_set_file | \
                        sort -t',' -k3 > "temp/${file_num}_finemap_filtered.txt"

                        # Check if there are any matching SNP ids if so then we join the files
                        if grep -qf temp/${file_num}_coloc_snps.txt <(cut -d',' -f3 "temp/${file_num}_finemap_filtered.txt"); then

                                # Perform join if matches are found

                                join -t ',' -1 3 -2 16 -o auto "temp/${file_num}_finemap_filtered.txt" \
                                "temp/catalogue_${file_num}_coloc_results.txt" >> "${output_array[$i]}"

			fi
                fi
	
	done

done

echo "Done!"
