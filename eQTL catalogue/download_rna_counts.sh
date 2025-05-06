#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

curl -O https://zenodo.org/api/records/4678936/files-archive
