#!/bin/bash
#SBATCH --mem=100G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

module load R/4.4.1

echo "Extracting gene ranks for eQTL catalogue"

#/nfs/sw/R/R-4.3.1/bin/Rscript 06.extract_eQTL_catalogue_sumstats.R

echo "Extracting gene ranks for GTEx"

#Rscript 06.gtex_gene_distances.R

echo "Extracting gene ranks for onek1k"

#Rscript 06.Onek1k_gene_distances.R

echo "Done"

Rscript 06.MAGE_gene_distances.R
