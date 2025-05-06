#!/bin/bash

array=("MON")

#array=("BAS" "EOS" "HCT" "HGB" "LYM" "MCHC" "MCH" "MCV" "MON" "MPV" "NEU" "PLT" "RBC" "RDW" "WBC")

module purge
module load R/4.3.1

echo 'Running SuSiE COLOC Rscript for each element in array'

# Loop through each element in the array and submit a job
for i in "${array[@]}"
do
  echo "Submitting job for ${i}"
  sbatch --mem=50G --cpus-per-task=10 --export=ALL,i="${i}" --wrap="Rscript multi_ancestry_coloc_beta.R ${i}"
done
