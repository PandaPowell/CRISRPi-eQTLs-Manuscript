#!/bin/bash

# This script merges all colocalised results by chromosome into a cell type file
# Then merges all cell type files into one overall file

mkdir -p coloc_results_combined


  	eqtl="GTEx"
	
	echo "SNP,hit2,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,chr.x,pos.x,a0.x,a1.x,minor_allele,minor_AF,low_confidence_variant,n_complete_samples,AC,ytx,beta.x,se,tstat,pval,Chr,Pos,_NUM_ID_.ss.x,rsid.x,_NUM_ID_.x,varbeta.x,z.x,chr.y,pos.y,a0.y,a1.y,snp_hg38.x,lower_hg38,upper_hg38,rsid.ss,lower_hg19,upper_hg19,phenotype_id,variant_id,tss_distance,af,ma_samples,ma_count,pval_nominal,beta.y,slope_se,_NUM_ID_.ss.y,rsid.y,_NUM_ID_.y,varbeta.y,z.y,trait,region,method,gwas" > coloc_results_combined/${eqtl}_coloc_results.txt

	for i in {1..22}; do
	
		chr="chr${i}"
	
		for file_num in $(seq -w 30000 10 30300); do
		
    		formatted_file="../results/Coloc_results_V2/${eqtl}/${file_num}/${chr}_coloc_results.txt"
		
    		if [ -f "$formatted_file" ]; then
    			
			
    			tail -n +2 $formatted_file >> coloc_results_combined/${eqtl}_coloc_results.txt
    	
		fi
		done
		
	done
	
	total=$(wc -l coloc_results_combined/${eqtl}_coloc_results.txt)
	
	echo "Total number of lines across all directories: $total"
