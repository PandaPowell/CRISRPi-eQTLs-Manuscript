#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=20G

# Define an array of strings
declare -a cell_types=("NK_cells" "B_cells" "CD4_T_cells" "CD8_T_cells" "DC_mean" "Mono_cells" "Other_cells" "Other_T_cells")

for cell_type in "${cell_types[@]}"
do

echo "Processing $cell_type"

base_dir="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/results/Coloc_results_V2/${cell_type}"


for file_num in $(seq -w 30000 10 30300); do
    dir="${base_dir}/${file_num}/"
    if [ -d "$dir" ]; then

        # Merge chr files for each GWAS

        echo $dir

	rm "${base_dir}/${file_num}_combined_coloc_results.txt"

	> "${base_dir}/${file_num}_combined_coloc_results.txt"

	for i in {1..22}
	do
		file="${dir}chr${i}_coloc_results.txt"
        	tail -n +2 "$file" >> "${base_dir}/${file_num}_combined_coloc_results.txt"
	done

    else
        echo "Directory does not exist: $dir"
    fi
done

# Merge each GWAS coloc results into one file by cell type

rm "${base_dir}/${cell_type}_combined_coloc_results.txt"

echo "SNP,hit2,PP.H0.abf,PP.H1.abf,PP.H2.abf,PP.H3.abf,PP.H4.abf,nsnps,gwas_credible_set,eqtl_credible_set,SNP.PP.H4,chr.x,pos.x,a0.x,a1.x,minor_allele,minor_AF,low_confidence_variant,n_complete_samples,AC,ytx,beta.x,se,tstat,pval,Chr,Pos,_NUM_ID_.ss.x,rsid.x,_NUM_ID_.x,varbeta.x,z.x,chr.y,pos.y,a0.y,a1.y,phenotype_id,variant_id,start_distance,af,ma_samples,ma_count,pval_nominal,beta.y,slope_se,_NUM_ID_.ss.y,rsid.y,_NUM_ID_.y,varbeta.y,z.y,trait,region,method,gwas" > "${base_dir}/${cell_type}_combined_coloc_results.txt"

for file_num in $(seq -w 30000 10 30300)
do

cat "${base_dir}/${file_num}_combined_coloc_results.txt" >> "${base_dir}/${cell_type}_combined_coloc_results.txt"

done

done
