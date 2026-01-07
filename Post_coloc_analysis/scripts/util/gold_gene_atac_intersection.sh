#!/bin/bash

outdir="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/ATAC_overlap"
mkdir -p "$outdir"

rm "$outdir/gold_genes_ATAC_peaks.bed"

for file in /gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/ATAC-seq/peaks/*; do
  filename=$(basename "$file")

  # skip files
  if [[ "$filename" == "README" || "$filename" == "other" ]]; then
    continue
  fi

  echo "Processing $file"

  # strip suffix to get cell name
  cell=${filename%%_ATAC_peaks.bed}

  # run bedtools and append cell name as extra column
  /nfs/sw/bedtools/bedtools-2.31.1/bin/bedtools intersect -loj \
    -a /gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/cres_w_grnas.bed \
    -b "$file" \
  | awk -v c="$cell" 'BEGIN{OFS="\t"} {print $0, c}' \
  >> "$outdir/gold_genes_ATAC_peaks.bed"

done
