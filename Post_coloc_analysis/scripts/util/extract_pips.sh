#!/bin/bash

for folder in $(seq 30000 10 30300); do
  for chr in {1..22}; do
    grep -Ff finemap_snp_intersect_grna.txt ../results/UKBB_SuSiE_finemap/${folder}/${folder}_chr${chr}_finemap_results_pip.0.01.txt
  done
done > extracted_pips.txt
