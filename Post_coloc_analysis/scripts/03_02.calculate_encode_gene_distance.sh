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

bedtools window -a ../CRISPR_data/sorted_GRC37.bed -b genes_tss.bed -w 1000000 > nearby_genes.bed

echo -e "chr\tgRNA_start\tgRNA_stop\ttarget_site\tTSS_start\tgene_stop\tgene\tdistance" > nearby_genes.fin.bed

awk -v OFS="\t" '{ start=$2; tss=$7; distance=(tss > start) ? tss-start : start-tss; print $1,$2,$3,$4,$7,$9,$8,distance }' nearby_genes.bed | \
awk -F"\t" '{gsub(/[";]/, "", $7); print}' OFS="\t" >> nearby_genes.fin.bed

#Filter to only genes expressed in K562 cells

awk '{print$9"\t"$2}' ../CRISPR_data/220308_STINGseq-CRE-Network-Genes.txt | sort -k1,1 > k562_genes.txt

awk 'NR==FNR {a[$1]; next} $7 in a' k562_genes.txt nearby_genes.fin.bed | sort -k7,7 > temp.bed

echo -e "gene_id\tchr\tgRNA_start\tgRNA_stop\ttarget_site\tTSS_start\tgene_stop\tdistance\tgene" > encode_k562_distances.bed

join -t $'\t' -1 7 -2 1 -o auto temp.bed k562_genes.txt | uniq >> encode_k562_distances.bed

echo "Done!"
