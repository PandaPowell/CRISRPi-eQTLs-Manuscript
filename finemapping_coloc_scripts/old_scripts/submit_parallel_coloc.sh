#!bin/bash

REGION_FILE=$1
GWAS_FILE=$2
GWAS_NAME="${GWAS_FILE##*/}"
GWAS_NAME2="${GWAS_NAME%.*}"
GWAS_NAME3="${GWAS_NAME2%_*}"

echo "FINEMAPPING GWAS" $GWAS_NAME2

echo "Spliting regions in region file by chromosome"

mkdir lead_loci
mkdir ./Reports/coloc_out

mkdir /gpfs/commons/groups/lappalainen_lab/sghatan/Coloc_results/${GWAS_NAME3}

for i in {1..22}
do
        awk -v chr=$i '{if ($1==chr) print $0}' $REGION_FILE > lead_loci/chr${i}.regions
done

for i in {1..22}
do
        sbatch submit_susie_coloc.sh $2 lead_loci/chr${i}.regions
done
