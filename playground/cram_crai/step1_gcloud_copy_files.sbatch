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

# -------------------------------
# Settings
# -------------------------------
CRAM_DIR=$SOL/cram
mkdir -p "$CRAM_DIR"

# -------------------------------
# Function: copy_crams_from_bucket
# -------------------------------
copy_crams_from_bucket () {
    local BUCKET_PATH=$1

    echo "Scanning bucket: $BUCKET_PATH"

    # List all cram files recursively
    gcloud storage ls -r "${BUCKET_PATH}" | grep '\.cram$' | while read -r CRAM_FILE; do
        FILE_NAME=$(basename "$CRAM_FILE")

        if [ ! -f "$CRAM_DIR/$FILE_NAME" ]; then
            echo "Copying $CRAM_FILE → $CRAM_DIR/"
            gcloud storage cp "$CRAM_FILE" "$CRAM_DIR/"
        else
            echo "Skipping $FILE_NAME (already exists)."
        fi
    done
}

# -------------------------------
# Example call
# -------------------------------
# Usage: sbatch this_script.sh gs://bucket-name/path/
if [ $# -eq 0 ]; then
    echo "Error: No Google bucket path provided."
    echo "Usage: sbatch $0 gs://your-bucket/path/"
    exit 1
fi

copy_crams_from_bucket "$1"

# GOOGLE_BUCKET='gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/'