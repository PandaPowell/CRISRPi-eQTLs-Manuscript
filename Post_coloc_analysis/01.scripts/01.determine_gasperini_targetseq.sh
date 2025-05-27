#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module purge
module load R/4.3.1
module load blast/2.10.0

cd "/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis"

mkdir -p gasperini_spaceseq

# Download enhancer positions and gRNA spacers
Rscript gasperini_spaceseq/download_format_supp.R

# Create blast database from reference genome
#makeblastdb -in $home/genome/GRCh37.p13.fasta/GCF_000001405.25/GCF_000001405.25_GRCh37.p13_genomic.fna \
#-title GRCh37.p13 -dbtype nucl -parse_seqids \
#-out $home/genome/GRCh37.p13.fasta/GCF_000001405.25/GRCh37.p13

# Create fasta file from spacer sequences
> "gasperini_spaceseq/gRNA_targetsite_spacers.fa"

tail -n +2 gasperini_spaceseq/gRNA_targetsite_spacers.txt | \
while read line
do
    CHR=$(echo $line | awk '{print $3}' | sed 's/chr//g')
    START=$(echo $line | awk '{print $4}')
    STOP=$(echo $line | awk '{print $5}')
    RANGE="$CHR:$START-$STOP"
    FASTA=$(echo $line | awk '{print $1}')
    TARGET=$(echo $line | awk '{print$8}')

    echo -e ">${TARGET}\n${FASTA}" >> "gasperini_spaceseq/gRNA_targetsite_spacers.fa"

done

echo "Fasta file created aligning to GRCh37.p13 using blast"

# Use blast to align 
blastn -db $home/genome/GRCh37.p13.fasta/GCF_000001405.25/GRCh37.p13 \
       -query gasperini_spaceseq/gRNA_targetsite_spacers.fa \
       -out gasperini_spaceseq/blast_results.out \
       -word_size 20 \
       -evalue 1 \
       -reward 1 \
       -penalty -3 \
       -gapopen 5 \
       -gapextend 2 \
       -dust no \
       -outfmt 6 \
       -num_threads 8 \
       -task blastn

echo "Running R script to clean blast output and make bed file for intersections"

Rscript gasperini_spaceseq/determine_spacer_seq.R

echo "Done!"
