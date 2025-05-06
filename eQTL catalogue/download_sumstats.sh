#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

DATASET_FILE=$1
REMOTE="http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/"

# Make sure the DATASET_FILE variable is correctly expanded
while read line; do

        STUDY_ID=$(echo $line | awk -F, '{print $1}')
        DATA_ID=$(echo $line | awk -F, '{print $2}')
	OUT_DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/${STUDY_ID}/${DATA_ID}/"

	rm "${OUT_DIR}${DATA_ID}.all.tsv.gz"

        echo "Downloading ${DATA_ID}..."

        # Use the correct variable (DATA_ID) instead of DATASET_ID
        wget -P ${OUT_DIR} "${REMOTE}${STUDY_ID}/${DATA_ID}/${DATA_ID}.all.tsv.gz"

        echo "Done!"

done < "$DATASET_FILE"
