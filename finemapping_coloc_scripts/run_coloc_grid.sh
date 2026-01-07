#!/bin/bash
#SBATCH --job-name=susie_coloc_R
#SBATCH --mem=30G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=7-00:00:00
#SBATCH -o ./Reports/coloc_out/logs/Interval_coloc-%A_%a.out

set -euo pipefail

# Inputs from launcher via --export (with sensible fallbacks for ad-hoc runs)
GWAS_LIST="${GWAS_LIST:-gwas_list.txt}"
N_CHR="${N_CHR:-22}"

# Map 1D array index -> (gwas_idx, chr)
task=${SLURM_ARRAY_TASK_ID}
if [[ -z "${task:-}" ]]; then
  echo "SLURM_ARRAY_TASK_ID not set." >&2
  exit 1
fi

gwas_idx=$(( (task - 1) / N_CHR + 1 ))
chr=$(( (task - 1) % N_CHR + 1 ))

# Fetch the gwas file for this index
if ! gwas_file=$(sed -n "${gwas_idx}p" "${GWAS_LIST}"); then
  echo "Failed to read GWAS index ${gwas_idx} from ${GWAS_LIST}" >&2
  exit 1
fi
if [[ -z "${gwas_file}" || ! -f "${gwas_file}" ]]; then
  echo "GWAS file missing for index ${gwas_idx}: ${gwas_file}" >&2
  exit 1
fi

module purge
module load R/4.4.1

echo "Running SuSiE COLOC Rscript for chr ${chr} and GWAS: ${gwas_file}"

NAME="Interval"
EQTL="/gpfs/commons/datasets/controlled/INTERVAL/public_sumstats/cis/INTERVAL_eQTL_summary_statistics/INTERVAL_eQTL_nominal_chr${chr}.tsv"

# Output dirs
mkdir -p "./Reports/coloc_out/${NAME}" "./Reports/coloc_out/logs"

# Optional: fail early if eQTL file is missing
if [[ ! -f "${EQTL}" ]]; then
  echo "Missing eQTL file: ${EQTL}" >&2
  exit 1
fi

# Call your R script:  NAME  chr  gwas_formatted  EQTL
Rscript 06.Interval_coloc.R "${NAME}" "${chr}" "${gwas_file}" "${EQTL}"

# ./launch_coloc_grid.sh
