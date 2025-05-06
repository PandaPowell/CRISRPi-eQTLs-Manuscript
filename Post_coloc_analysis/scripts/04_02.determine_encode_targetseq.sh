#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module purge
module load R/4.3.1
module load blast/2.10.0
module load liftover

mkdir -p encode_spaceseq
mkdir -p CRISPR_data

# Obtain harmonised CRISPRi data from encode
#wget -O CRISPR_data/ENCODE_harmonised_CRISPRi_data_GRC38.tsv https://www.encodeproject.org/files/ENCFF968BZL/@@download/ENCFF968BZL.tsv

# Create blast database from reference genome
makeblastdb -in $home/genome/GRCh38/GCF_000001405.26_GRCh38_genomic.fna \
-title GRCh38 -dbtype nucl -parse_seqids \
-out $home/genome/GRCh38/GRCh38

# Create fasta file from spacer sequences
> "encode_spaceseq/gRNA_targetsite_spacers.fa"

tail -n +2 CRISPR_data/ENCODE_harmonised_CRISPRi_data_GRC38.tsv | \
while read line
do
    CHR=$(echo $line | awk '{print $1}' | sed 's/chr//g')
    START=$(echo $line | awk '{print $2}')
    STOP=$(echo $line | awk '{print $3}')
    RANGE="$CHR:$START-$STOP"
    FASTA=$(echo $line | awk '{split($17,a,";"); print a[1]}')
    TARGET=$(echo $line | awk '{print$7}')

    echo -e ">${TARGET}\n${FASTA}" >> "encode_spaceseq/gRNA_targetsite_spacers.fa"

done

echo "Fasta file created aligning to GRCh38 using blast"

# Use blast to align
blastn -db $home/genome/GRCh38/GRCh38 \
       -query encode_spaceseq/gRNA_targetsite_spacers.fa \
       -out encode_spaceseq/blast_results.out \
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

Rscript encode_spaceseq/determine_spacer_seq.R

echo "Alignment Done!"

echo "Now lifting over to GRC37"

cd CRISPR_data

# Liftover to GRC37
liftover ../encode_spaceseq/encode_gRNA_target_sites_GRCh38.bed \
$home/stingseq_eqtl_overlap/data/liftover_chain_files/hg38ToHg19.over.chain.gz \
ENCODE_harmonised_CRISPRi_data_GRC37.bed unmapped_regions

cat unmapped_regions

# Sort files before joining
sort -k4,4 ENCODE_harmonised_CRISPRi_data_GRC37.bed > sorted_GRC37.bed
sed 's/ \+/\t/g' ENCODE_harmonised_CRISPRi_data_GRC38.tsv | \
sort -k7,7 > sorted_GRC38.tsv

# Make header
echo -e "target_site\tchrom_GRC37\tchromStart_GRC37\tchromEnd_GRC37\tchrom_GRC38\tchromStart_GRC38\tchromEnd_GRC38\tname\tEffectsize\tstrandPerturbationTarget\tchrTSS\tstartTSS\tendTSS\tstrandGene\tEffectSize95ConfidenceIntervalLow\tEffectSize95ConfidenceIntervalHigh\tmeasuredGeneSymbol\tmeasuredEnsemblID\tguideSpacerSeq\tguideSeq\tSignificant\tpValue\tpValueAdjusted\tPowerAtEffectSize25\tPowerAtEffectSize10\tPowerAtEffectSize15\tPowerAtEffectSize20\tPowerAtEffectSize50\tValidConnection\tNotes\tReference#" > ENCODE_harmonised_CRISPRi_data_GRC37.tsv

# Join the sorted files based on the 4th column
join -t $'\t' -1 4 -2 7 -o auto sorted_GRC37.bed sorted_GRC38.tsv >> ENCODE_harmonised_CRISPRi_data_GRC37.tsv

grep -v Gasperini ENCODE_harmonised_CRISPRi_data_GRC37.tsv > NoGasperini_crispri_data.tsv

echo "Download, liftover and formatting Done!"
