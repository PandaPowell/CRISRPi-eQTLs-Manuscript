#!/bin/bash

module load bedtools/2.31.0-GCC-12.3.0

# Format cres overlapping gRNAs to bed file
tail -n +2 cres_with_grnas.txt | awk -F, '{print "chr"$11"\t"$13-1"\t"$13"\t"$12}' | sort -u | sort -k1,1 -k2,2n > cres_w_grnas.bed

# Filter K562 ATAC-seq peaks to only those that overlap a gRNA
bedtools closest -d -a cres_w_grnas.bed -b ../data/ATAC-seq/peaks/K562_ATAC_peaks.bed | awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$4"\t"$9}' \
> ATAC_overlap/K562_ATAC_peaks_w_cres.bed

echo -e "chr\tk562_pos_lower\tk562_pos_upper\tk562_peak\tgRNA_target\tdistance\tchr.x\tcell_pos_lower\tcell_pos_upper\tcell_peak" \
> ATAC_overlap/intersecting_ATAC_peaks_w_cres.bed

for file in ../data/ATAC-seq/peaks/*; do
  filename=$(basename "$file")
  
  if [ "$filename" = "K562_ATAC_peaks.bed" ]; then
    continue
  fi
  
  echo "Processing $file"
  # Filter K562 peaks with peaks from other cell types
  bedtools intersect -loj -a ATAC_overlap/K562_ATAC_peaks_w_cres.bed -b "$file" \
  >> ATAC_overlap/intersecting_ATAC_peaks_w_cres.bed

done

echo -e "chr\tk562_pos_lower\tk562_pos_upper\tk562_peak\tchr.x\tcell_pos_lower\tcell_pos_upper\tcell_peak" \
> ATAC_overlap/intersecting_ATAC_peaks.bed

for file in ../data/ATAC-seq/peaks/*; do
  filename=$(basename "$file")

  if [ "$filename" = "K562_ATAC_peaks.bed" ]; then
    continue
  fi

  echo "Processing $file"
  
  # Filter K562 peaks with peaks from other cell types
  bedtools intersect -loj -a ../data/ATAC-seq/peaks/K562_ATAC_peaks.bed -b "$file" \
  >> ATAC_overlap/intersecting_ATAC_peaks.bed

done
