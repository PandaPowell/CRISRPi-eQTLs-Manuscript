#!/bin/bash

cell_types=("NK_cells" "B_cells" "CD4_T_cells" "CD8_T_cells" "DC_mean" "Other_T_cells" "Other_cells")

path=("NK_mean_mx/Cis_GW_allpairs/OneK1K_NK.cis_qtl_pairs" "B_mean_mx/Cis_GW_allpairs/OneK1K_B.cis_qtl_pairs" "CD4_T_mean_mx/Cis_GW_allpairs/OneK1K_CD4_T.cis_qtl_pairs" \
"CD8_T_mean_mx/Cis_GW_allpairs/OneK1K_CD8_T.cis_qtl_pairs" "DC_mean_mx/Cis_GW_allpairs/OneK1K_DC.cis_qtl_pairs" \
"other_T_mean_mx/Cis_GW_allpairs/OneK1K_other_T.cis_qtl_pairs" "other_mean_mx/Cis_GW_allpairs/OneK1K_other.cis_qtl_pairs")

for j in $(seq 0 6)
do

NAME="${cell_types[$j]}"
EQTL="${path[$j]}"

echo "${NAME} : ${EQTL}"

sbatch 03.submit_susie_coloc.sh ${NAME} ${EQTL}

done
