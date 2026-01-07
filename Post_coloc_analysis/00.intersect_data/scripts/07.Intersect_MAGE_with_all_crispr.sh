#!/bin/bash

# set dir
cd "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

mkdir -p temp

echo "GWAS_variant,sentenial_variant,chr,finemap_snp_intersect_grna,SS_coord,eQTL_variant,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,SNP_coords.x,rs_number,reference_allele,other_allele,eaf,beta,se,beta_95L,beta_95U,z.x,p-value,_-log10_p-value,q_statistic,q_p-value,i2,n_studies,n_samples,effects,chr.1,pos,varbeta.x,SNP_coords.y,snp_hg38,chr.2,lower_hg38,upper_hg38,rsid,lower_hg19,upper_hg19,variantChrom,variantPosition,variant_kgpID,variant_rsID,ensemblID,geneSymbol,tss_distance,ma_samples,ma_count,pval_nominal,slope,slope_se,ref_h,alt_h,raf,varbeta.y,z.y,trait,region,method,gwas" \
>  "stingseq_eqtl_overlap/MAGE_coloc_stingseq.txt"

echo "GWAS_variant,sentenial_variant,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,eQTL_variant,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,SNP_coords.x,rs_number,reference_allele,other_allele,eaf,beta,se,beta_95L,beta_95U,z.x,p-value,_-log10_p-value,q_statistic,q_p-value,i2,n_studies,n_samples,effects,chr.1,pos,varbeta.x,SNP_coords.y,snp_hg38,chr.2,lower_hg38,upper_hg38,rsid,lower_hg19,upper_hg19,variantChrom,variantPosition,variant_kgpID,variant_rsID,ensemblID,geneSymbol,tss_distance,ma_samples,ma_count,pval_nominal,slope,slope_se,ref_h,alt_h,raf,varbeta.y,z.y,trait,region,method,gwas" \
> "gasperini_eqtl_overlap/MAGE_coloc_gasperini.txt"

echo "GWAS_variant,sentenial_variant,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,eQTL_variant,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,SNP_coords.x,rs_number,reference_allele,other_allele,eaf,beta,se,beta_95L,beta_95U,z.x,p-value,_-log10_p-value,q_statistic,q_p-value,i2,n_studies,n_samples,effects,chr.1,pos,varbeta.x,SNP_coords.y,snp_hg38,chr.2,lower_hg38,upper_hg38,rsid,lower_hg19,upper_hg19,variantChrom,variantPosition,variant_kgpID,variant_rsID,ensemblID,geneSymbol,tss_distance,ma_samples,ma_count,pval_nominal,slope,slope_se,ref_h,alt_h,raf,varbeta.y,z.y,trait,region,method,gwas" \
> "encode_eqtl_overlap/MAGE_coloc_crispri.txt"

array=("BAS" "EOS" "HCT" "HGB" "LYM" "MCHC" "MCH" "MCV" "MON" "MPV" "NEU" "PLT" "RBC" "RDW" "WBC")

for gwas in "${array[@]}"; do

	echo "Running for ${gwas}"

	coloc_file="$home/stingseq_eqtl_overlap/results/Coloc_results_V2/MAGE/MAGE_${gwas}_coloc_res.txt"
	credset_files=("$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data/sting_seq_credible_sets/BCX/${gwas}_credset.txt" "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data/gasperini_credible_sets/BCX/${gwas}_credset.txt" "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data/encode_credible_sets/BCX/${gwas}_credset.txt")
	output_array=("stingseq_eqtl_overlap/MAGE_coloc_stingseq.txt" "gasperini_eqtl_overlap/MAGE_coloc_gasperini.txt" "encode_eqtl_overlap/MAGE_coloc_crispri.txt")

	# Extract unique SNP ids from column of coloc file
        awk -F, '{sub(/^chr/, "", $1); print $1}' ${coloc_file} | \
        sort -u > temp/${gwas}_coloc_snps.txt

        # Format and sort coloc file
        awk -v OFS="," -F, '{sub(/^chr/, "", $1); print}' ${coloc_file} | \
        sort -t',' -k1 > "temp/mage_${gwas}_coloc_results.txt"

	for i in "${!credset_files[@]}"; do

		echo "Running for ${output_array[$i]}"

		if [ "${credset_files[$i]}" = "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data/sting_seq_credible_sets/BCX/${gwas}_credset.txt" ]; then

                        # Filter stingseq credible set file to include only lines where SNP from coloc file is present
                        awk -v OFS=',' 'NR==FNR {vals[$1]; next} $3 in vals {print $1, $2, $3, $4, $5}' "temp/${gwas}_coloc_snps.txt" \
                        "${credset_files[$i]}" | \
                        sort -t',' -k3 > "temp/${gwas}_finemap_filtered.txt"

                        # Check if there are any matching SNP ids if so then we join the files
                        if grep -qf temp/${gwas}_coloc_snps.txt <(cut -d',' -f3 "temp/${gwas}_finemap_filtered.txt"); then

                                # Perform join if matches are found

                                join -t ',' -1 3 -2 1 -o auto "temp/${gwas}_finemap_filtered.txt" \
                                "temp/mage_${gwas}_coloc_results.txt" >> "${output_array[$i]}"
                        
			fi
		
		else
	
                        # Filter stingseq credible set file to include only lines where SNP from coloc file is present
                       	awk -v OFS=',' 'NR==FNR {vals[$1]; next} $3 in vals {print $1, $2, $3, $4, $5, $6, $7}' "temp/${gwas}_coloc_snps.txt" \
                       	"${credset_files[$i]}" | \
                        sort -t',' -k3 > "temp/${gwas}_finemap_filtered.txt"

                        # Check if there are any matching SNP ids if so then we join the files
                       	if grep -qf temp/${gwas}_coloc_snps.txt <(cut -d',' -f3 "temp/${gwas}_finemap_filtered.txt"); then

                               	# Perform join if matches are found

                               	join -t ',' -1 3 -2 1 -o auto "temp/${gwas}_finemap_filtered.txt" \
                               	"temp/mage_${gwas}_coloc_results.txt" >> "${output_array[$i]}"

                        fi
                fi

	done

done

echo "Done!"
