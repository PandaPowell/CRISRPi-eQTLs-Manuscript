#!/bin/bash
#SBATCH --mem=20G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8

wget -O GSE120861_at_scale_screen.exprs.mtx.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE120861&format=file&file=GSE120861_at_scale_screen.exprs.mtx.gz"
