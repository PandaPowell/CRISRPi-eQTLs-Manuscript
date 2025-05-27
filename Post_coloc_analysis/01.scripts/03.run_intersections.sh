#!/bin/bash

cd /gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/Post_coloc_analysis/00.intersect_data/scripts/

# Make the scripts executable (if not already executable)
chmod +x 01.Filter_finemap_pip.sh \
	 01.Intersect_ENCODE_w_finemap.sh \
         01.Intersect_BCX.sh \
         01.Intersect_sting_seq_w_finemap.sh \
         02.Intersect_gasperini_w_finemap_closest.sh \
         03.Intersect_onek1k_gtex_with_stingseq.sh \
         04.Intersect_onek1k_gtex_with_gasperini.sh \
         05.Intersect_onek1k_gtex_with_encode.sh \
         06.Intersect_eqtl_cat_with_all_crispr.sh \
         07.Intersect_MAGE_with_all_crispr.sh

# Run scripts sequentially and stop if any script fails
set -e  # Exit immediately if a command exits with a non-zero status

echo "Running 01.Filter_finemap_pip.sh"
./01.Filter_finemap_pip.sh

echo "Running 00.Intersect_ENCODE_w_finemap.sh..."
./01.Intersect_ENCODE_w_finemap.sh

echo "Running 01.Intersect_BCX.sh..."
./01.Intersect_BCX.sh

echo "Running 01.Intersect_sting_seq_w_finemap.sh..."
./01.Intersect_sting_seq_w_finemap.sh

echo "Running 02.Intersect_gasperini_w_finemap_closest.sh..."
./02.Intersect_gasperini_w_finemap_closest.sh

echo "Running 03.Intersect_onek1k_gtex_with_stingseq.sh..."
./03.Intersect_onek1k_gtex_with_stingseq.sh

echo "Running 04.Intersect_onek1k_gtex_with_gasperini.sh..."
./04.Intersect_onek1k_gtex_with_gasperini.sh

echo "Running 05.Intersect_onek1k_gtex_with_encode.sh..."
./05.Intersect_onek1k_gtex_with_encode.sh

echo "Running 06.Intersect_eqtl_cat_with_all_crispr.sh..."
./06.Intersect_eqtl_cat_with_all_crispr.sh

echo "Running 07.Intersect_MAGE_with_all_crispr.sh..."
./07.Intersect_MAGE_with_all_crispr.sh

echo "All scripts executed successfully!"
