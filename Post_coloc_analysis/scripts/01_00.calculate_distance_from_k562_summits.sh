#!/bin/bash

module purge
module load bedtools

mkdir -p data
cd data

# Download ensemble data
#wget https://ftp.ensembl.org/pub/grch37/release-112/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz

#gunzip Homo_sapiens.GRCh37.87.gtf.gz

# Extract TSS positions
awk '$3 == "transcript" {print "chr"$1 "\t" ($4-1) "\t" $4 "\t" $10 "\t" $5 "\t" $7}' Homo_sapiens.GRCh37.87.gtf | \
sort -k1,1 -k2,2n | awk '!seen[$4]++' > genes_tss.bed

bedtools window -w 1000000 \
-a /gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/ATAC-seq/summits/K562_ATAC_summits.bed \
-b genes_tss.bed | \
awk '{print$1,$3,$4,$7,$8}' OFS='\t' | awk -F"\t" '{gsub(/[";]/, "", $5); print}' OFS="\t" > gene_dis_summits.bed

#Filter to only genes expressed in K562 cells

awk '{print$9"\t"$2}' ../CRISPR_data/220308_STINGseq-CRE-Network-Genes.txt | sort -k1,1 > k562_genes.txt

echo -e "chr\tsummit_pos\tsummit\tTSS_pos\tgene" > k562_gene_dis_summits.bed

awk 'NR==FNR {a[$1]; next} $5 in a' k562_genes.txt gene_dis_summits.bed >> k562_gene_dis_summits.bed

echo "Done!"
