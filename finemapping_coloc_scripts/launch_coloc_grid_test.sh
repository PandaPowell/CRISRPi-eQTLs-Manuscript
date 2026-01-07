#!/usr/bin/env bash
set -euo pipefail

# Where your formatted GWAS files live
BASE_DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/data/UKBB_sumstats"

# Range of candidate GWAS IDs to check (inclusive)
START_ID=30000
END_ID=30010

# Build a clean list of actual files
GWAS_LIST="gwas_list.txt"
: > "${GWAS_LIST}"
for id in $(seq -w ${START_ID} ${END_ID}); do
  f="${BASE_DIR}/${id}_formatted.tsv"
  [[ -f "$f" ]] && echo "$f" >> "${GWAS_LIST}"
done

N_GWAS=$(wc -l < "${GWAS_LIST}")
if [[ "${N_GWAS}" -eq 0 ]]; then
  echo "No GWAS formatted files found in ${BASE_DIR} for ${START_ID}-${END_ID}" >&2
  exit 1
fi

N_CHR=2
TOTAL=$(( N_GWAS * N_CHR ))

echo "Submitting ${TOTAL} array tasks for ${N_GWAS} GWAS × ${N_CHR} chromosomes."

# Concurrency cap (optional): change %100 to taste or remove it
sbatch --array=1-${TOTAL}%100 --export=ALL,GWAS_LIST="${GWAS_LIST}",N_CHR="${N_CHR}" run_coloc_grid_test.sh
