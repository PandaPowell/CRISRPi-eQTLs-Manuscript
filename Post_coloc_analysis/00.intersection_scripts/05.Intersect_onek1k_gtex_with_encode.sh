#!/bin/bash

# set dir
cd "$home/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data"

mkdir -p temp

mkdir -p encode_eqtl_overlap

cell_types=("NK_cells" "B_cells" "CD4_T_cells" "CD8_T_cells" "DC_mean" "Other_T_cells" "Mono_cells" "Other_cells")

for j in $(seq 0 7)
do
	eqtl="${cell_types[$j]}"

	echo "GWAS_variant,region,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,hit2,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,chr.x,pos.x,a0.x,a1.x,minor_allele,minor_AF,low_confidence_variant,n_complete_samples,AC,ytx,beta.x,se,tstat,pval,Chr,Pos,_NUM_ID_.ss.x,rsid.x,_NUM_ID_.x,varbeta.x,z.x,chr.y,pos.y,a0.y,a1.y,phenotype_id,variant_id,start_distance,af,ma_samples,ma_count,pval_nominal,beta.y,slope_se,_NUM_ID_.ss.y,rsid.y,_NUM_ID_.y,varbeta.y,z.y,trait,region.x,method,gwas_name" > "encode_eqtl_overlap/${eqtl}_coloc_crispri.txt"

	for i in {1..22}
	do
		
		chr="chr${i}"
	
		for file_num in $(seq -w 30000 10 30300); do
		
        		formatted_file="../../results/Coloc_results_V2/${eqtl}/${file_num}/${chr}_coloc_results.txt"
		
        		if [ -f "$formatted_file" ]; then
		
			# Extract unique GWAS snps
			cut -d',' -f1 $formatted_file | \
			sort -u > temp/${file_num}_${chr}_coloc_snps.txt
		
			# Filter finemapped snps that intersect gasperini targets to those that colocalise an eqtl
			awk -v OFS=',' 'NR==FNR {vals[$1]; next} $3 in vals {print $1, $2, $3, $4, $5, $6, $7}' "temp/${file_num}_${chr}_coloc_snps.txt" \
			"encode_credible_sets/${file_num}/${file_num}_${chr}_stingseq_credset.txt" | \
			sort -t',' -k3 > "temp/${file_num}_${chr}_finemap_filtered.txt"
		
			# Sort GWAS coloc file
			tail -n +2 $formatted_file | sort -t',' -k1 > "temp/${eqtl}_${file_num}_${chr}_coloc_results.txt"
		
			# Check if there are any matching snps, if so then join the files
			if grep -qf temp/${file_num}_${chr}_coloc_snps.txt <(cut -d',' -f3 "temp/${file_num}_${chr}_finemap_filtered.txt"); then
	   	 		# Perform join if matches are found

	    			join -t ',' -1 3 -2 1 -o auto "temp/${file_num}_${chr}_finemap_filtered.txt" \
	    			"temp/${eqtl}_${file_num}_${chr}_coloc_results.txt" >> "encode_eqtl_overlap/${eqtl}_coloc_crispri.txt"

			fi
        		fi
		done
	done	
done

echo "GWAS_variant,region_cs,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,eQTL_variant,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,chr.x,pos.x,a0.x,a1.x,minor_allele,minor_AF,low_confidence_variant,n_complete_samples,AC,ytx,beta.x,se,tstat,pval,Chr,Pos,_NUM_ID_.ss.x,rsid.x,_NUM_ID_.x,varbeta.x,z.x,chr.y,pos.y,a0.y,a1.y,phenotype_id,variant_id,start_distance,af,ma_samples,ma_count,pval_nominal,beta.y,slope_se,_NUM_ID_.ss.y,rsid.y,_NUM_ID_.y,varbeta.y,z.y,trait,region,method,gwas_name,eqtl" \
> "encode_eqtl_overlap/All_coloc_crispri_overlap.txt"

for j in $(seq 0 7)
do
  	eqtl="${cell_types[$j]}"

        echo "encode_eqtl_overlap/${eqtl}_coloc_crispri.txt"

        tail -n +2 "encode_eqtl_overlap/${eqtl}_coloc_crispri.txt" | \
	awk -v OFS=',' -v eqtl="$eqtl" '{print $0, eqtl}' >> "encode_eqtl_overlap/All_coloc_crispri_overlap.txt"

done

echo "Intersecting GTEx coloc"

eqtl="GTEx"

echo "GWAS_variant,region,chr,finemap_snp_intersect_grna,lower_grna_target_site,upper_grna_target_site,target_site,eQTL_variant,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,chr.x,pos.x,a0.x,a1.x,minor_allele,minor_AF,low_confidence_variant,n_complete_samples,AC,ytx,beta.x,se,tstat,pval,Chr,Pos,_NUM_ID_.ss.x,rsid.x,_NUM_ID_.x,varbeta.x,z.x,chr.y,pos.y,a0.y,a1.y,snp_hg38.x,lower_hg38,upper_hg38,rsid.ss,lower_hg19,upper_hg19,phenotype_id,variant_id,tss_distance,af,ma_samples,ma_count,pval_nominal,beta.y,slope_se,_NUM_ID_.ss.y,rsid.y,_NUM_ID_.y,varbeta.y,z.y,trait,region.x,method,gwas_name" > "encode_eqtl_overlap//${eqtl}_coloc_crispri.txt"

        for i in {1..22}
        do

                chr="chr${i}"

                for file_num in $(seq -w 30000 10 30300); do

                        formatted_file="../../results/Coloc_results_V2/${eqtl}/${file_num}/${chr}_coloc_results.txt"

                        if [ -f "$formatted_file" ]; then

                        # Extract unique GWAS snps
                        cut -d',' -f1 $formatted_file | \
                        sort -u > temp/${file_num}_${chr}_coloc_snps.txt

                        # Filter finemapped snps that intersect gasperini targets to those that colocalise an eqtl
                        awk -v OFS=',' 'NR==FNR {vals[$1]; next} $3 in vals {print $1, $2, $3, $4, $5, $6, $7}' "temp/${file_num}_${chr}_coloc_snps.txt" \
                        "encode_credible_sets/${file_num}/${file_num}_${chr}_stingseq_credset.txt" | \
                        sort -t',' -k3 > "temp/${file_num}_${chr}_finemap_filtered.txt"

                        # Sort GWAS coloc file
                        tail -n +2 $formatted_file | sort -t',' -k1 > "temp/${eqtl}_${file_num}_${chr}_coloc_results.txt"

                        # Check if there are any matching snps, if so then join the files
                        if grep -qf temp/${file_num}_${chr}_coloc_snps.txt <(cut -d',' -f3 "temp/${file_num}_${chr}_finemap_filtered.txt"); then
                        
                                # Perform join if matches are found
                                
                                join -t ',' -1 3 -2 1 -o auto "temp/${file_num}_${chr}_finemap_filtered.txt" \
                                "temp/${eqtl}_${file_num}_${chr}_coloc_results.txt" >> "encode_eqtl_overlap/${eqtl}_coloc_crispri.txt"

                        fi
                        fi
                done
        done


rm -r temp
