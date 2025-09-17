#!/bin/bash
#SBATCH --job-name=gcloud_copy
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=5G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

module purge
module load gcloud

CRAM_DIR=$SOL/cram 
GOOGLE_BUCKET='gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/'

# Create local directory if it doesn't exist
mkdir -p "$CRAM_DIR"

# Iterate through all .cram files in the bucket recursively
gcloud storage ls -r "${GOOGLE_BUCKET}" | grep '\.cram$' | while read -r CRAM_FILE; do
    # Get the base filename
    FILE_NAME=$(basename "$CRAM_FILE")

    # Only copy if file doesn't already exist locally
    if [ ! -f "$CRAM_DIR/$FILE_NAME" ]; then
        echo "Copying $CRAM_FILE to $CRAM_DIR/"
        gcloud storage cp "$CRAM_FILE" "$CRAM_DIR/"
    else
        echo "Skipping $FILE_NAME; already exists."
    fi
done
