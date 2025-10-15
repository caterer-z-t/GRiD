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
# USAGE
# -------------------------------
usage() {
    echo "Usage: sbatch $0 [-d CRAM_DIR] BUCKET_PATH"
    echo ""
    echo "  BUCKET_PATH   Google Cloud bucket path (required, e.g. gs://your-bucket/path/)"
    echo "  -d CRAM_DIR   Local directory to copy cram files into (default: \$SOL/cram)"
    echo "  -h            Show this help message"
    exit 1
}

# -------------------------------
# Defaults
# -------------------------------
CRAM_DIR=$SOL/cram

# -------------------------------
# Parse arguments
# -------------------------------
while getopts "d:h" opt; do
    case $opt in
        d) CRAM_DIR="$OPTARG" ;;
        h) usage ;;
        \?) usage ;;
    esac
done
shift $((OPTIND - 1))   # shift past parsed flags

# -------------------------------
# Validate arguments
# -------------------------------
if [ $# -eq 0 ]; then
    usage "Error: Missing BUCKET_PATH"
fi

BUCKET_PATH=$1

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

remove_empty_files() {
    find "$CRAM_DIR" -type f -empty -exec rm -f {} \;
}


# -------------------------------
# Run
# -------------------------------
copy_crams_from_bucket "$BUCKET_PATH"
remove_empty_files
