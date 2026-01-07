#!/bin/bash

module purge
module load bedtools

# Download ABC predictions
#wget ftp://ftp.broadinstitute.org/outgoing/lincRNA/ABC/AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt.gz

# Make array of cell type names in ABC file
readarray -t cell_types < <(cut -f$(head -n1 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt | awk -F'\t' '{print NF}') \
  AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt \
  | tail -n +2 \
  | sort -u \
  | grep -iE '(^|_)(B_cell|T_cell|T-cell|monocyte|macrophage|dendritic|natural_killer|NK|lymphocyte|CD[0-9]+|GM12878|K562|Jurkat|BJAB|U937|OCI|erythroblast|megakaryocyte|spleen|thymus)([^a-zA-Z]|$)')

array_length=${#cell_types[@]}

DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis"

# make intersections file header
echo -en "chr\tlower_grna_pos\tupper_grna_pos\ttarget_site\t" > "cres_w_grnas_ABC_intersections.txt"
head -n 1 AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt >> "cres_w_grnas_ABC_intersections.txt"

# Iterate over the array using the length
for (( i=0; i<array_length; i++ )); do
  
  # Split file into cell types
  name=${cell_types[$i]}
  echo "Element $i: ${name}"
  grep ${name} AllPredictions.AvgHiC.ABC0.015.minus150.ForABCPaperV3.txt > "${name}_ABC_predictions.txt"

  bedtools intersect -wb -a "cres_w_grnas.bed" -b "${name}_ABC_predictions.txt" >> "cres_w_grnas_ABC_intersections.txt"

done

echo "Done!"
