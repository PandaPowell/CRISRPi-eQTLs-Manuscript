module purge
module load liftover

# Obtain harmonised CRISPRi data from encode
#wget -O ENCODE_harmonised_CRISPRi_data_GRC38.tsv https://www.encodeproject.org/files/ENCFF968BZL/@@download/ENCFF968BZL.tsv

# Format file for liftover
tail -n +2 ENCODE_harmonised_CRISPRi_data_GRC38.tsv | awk '{print $1, $2, $3, $4}' > ENCODE_harmonised_CRISPRi_data_GRC38.bed

# Liftover to GRC37
liftover ENCODE_harmonised_CRISPRi_data_GRC38.bed \
$home/stingseq_eqtl_overlap/data/liftover_chain_files/hg38ToHg19.over.chain.gz \
ENCODE_harmonised_CRISPRi_data_GRC37.bed unmapped_regions

cat unmapped_regions

# Sort files before joining
sort -k4,4 ENCODE_harmonised_CRISPRi_data_GRC37.bed > sorted_GRC37.bed
sed 's/ \+/\t/g' ENCODE_harmonised_CRISPRi_data_GRC38.tsv | \
sort -k4,4 > sorted_GRC38.tsv

# Make header
echo -e "name\tchrom_GRC37\tchromStart_GRC37\tchromEnd_GRC37\tchrom_GRC38\tchromStart_GRC38\tchromEnd_GRC38\tEffectSize\tstrandPerturbationTarget\tPerturbationTargetID\tchrTSS\tstartTSS\tendTSS\tstrandGene\tEffectSize95ConfidenceIntervalLow\tEffectSize95ConfidenceIntervalHigh\tmeasuredGeneSymbol\tmeasuredEnsemblID\tguideSpacerSeq\tguideSeq\tSignificant\tpValue\tpValueAdjusted\tPowerAtEffectSize25\tPowerAtEffectSize10\tPowerAtEffectSize15\tPowerAtEffectSize20\tPowerAtEffectSize50\tValidConnection\tNotes\tReference#" > ENCODE_harmonised_CRISPRi_data_GRC37.tsv

# Join the sorted files based on the 4th column
join -t $'\t' -1 4 -2 4 -o auto sorted_GRC37.bed sorted_GRC38.tsv >> ENCODE_harmonised_CRISPRi_data_GRC37.tsv

grep -v Gasperini ENCODE_harmonised_CRISPRi_data_GRC37.tsv > NoGasperini_crispri_data.tsv

echo "Done!"
