#!/bin/bash

module load bedtools/2.31.0-GCC-12.3.0

# 1) Build cres_w_grnas.bed (unique by $12) and sort it (BED: chr, start, end)
tail -n +2 processed_data/cres_with_grnas.txt \
| awk -F',' 'BEGIN{OFS="\t"}
{
  k = $12
  if (!(k in seen)) {
    seen[k] = 1
    # BED is 0-based, half-open; make 1-bp interval at 1-based position $13
    print "chr"$11, $13-1, $13, k
  }
}' \
| LC_ALL=C sort -k1,1 -k2,2n -k3,3n \
> cres_w_grnas.bed

# Filter K562 ATAC-seq peaks to only those that overlap a gRNA
bedtools closest -d -a cres_w_grnas.bed -b ../data/ATAC-seq/peaks/K562_ATAC_peaks.bed | awk '{print $5"\t"$6"\t"$7"\t"$8"\t"$4"\t"$9}' \
> ATAC_overlap/K562_ATAC_peaks_w_cres.bed

echo -e "chr\tlower\tupper\tgRNA_target\tchr.x\tcell_pos_lower\tcell_pos_upper\tcell_peak" \
> ATAC_overlap/intersecting_ATAC_peaks_w_cres.bed

# 4) Loop over other peak files; skip selected names safely
for file in ../data/ATAC-seq/peaks/*; do
  filename=$(basename "$file")

  case "$filename" in
    K562_ATAC_peaks.bed|other|README)  # skip these
      continue
      ;;
  esac

  echo "Processing $file"

  # Ensure B is sorted; write to a temp in /tmp to avoid clobbering originals
  tmp_b=$(mktemp)
  LC_ALL=C sort -k1,1 -k2,2n -k3,3n "$file" > "$tmp_b"

  # Intersect (left outer join), append
  bedtools intersect -sorted -loj \
    -a cres_w_grnas.bed \
    -b "$tmp_b" \
  >> ATAC_overlap/intersecting_ATAC_peaks_w_cres.bed

  rm -f "$tmp_b"
done

echo -e "chr\tk562_pos_lower\tk562_pos_upper\tk562_peak\tchr.x\tcell_pos_lower\tcell_pos_upper\tcell_peak" \
> ATAC_overlap/intersecting_ATAC_peaks.bed

for file in ../data/ATAC-seq/peaks/*; do
  filename=$(basename "$file")

  case "$filename" in
    K562_ATAC_peaks.bed|other|README)  # skip these
      continue
      ;;
  esac

  echo "Processing $file"
  
  # Filter K562 peaks with peaks from other cell types
  bedtools intersect -loj -a ../data/ATAC-seq/peaks/K562_ATAC_peaks.bed -b "$file" \
  >> ATAC_overlap/intersecting_ATAC_peaks.bed

done
