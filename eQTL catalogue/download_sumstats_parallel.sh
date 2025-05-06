#!/bin/bash
#SBATCH --mem=20G
#SBATCH --cpus-per-task=8

DATASET_FILE=$1

# Function to handle downloading for each line
download_data() {
    
    line=$1
    STUDY_ID=$(echo "$line" | awk -F, '{print $1}')
    DATA_ID=$(echo "$line" | awk -F, '{print $2}')
    OUT_DIR="/gpfs/commons/groups/lappalainen_lab/sghatan/stingseq_eqtl_overlap/eQTL_catalogue/${STUDY_ID}/${DATA_ID}/"

    # Remove existing file if it exists
    rm -f "${OUT_DIR}${DATA_ID}.all.tsv.gz"

    echo "Downloading ${DATA_ID}..."

    # Download the file
    wget -P "${OUT_DIR}" "http://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/${STUDY_ID}/${DATA_ID}/${DATA_ID}.all.tsv.gz"

    echo "Done with ${STUDY_ID}${DATA_ID}!"
}

export -f download_data  # Export the function so parallel can use it

# Run the download function in parallel using all available CPUs
tail -n +2 "$DATASET_FILE" | parallel --keep-order -j $SLURM_CPUS_PER_TASK download_data {}
