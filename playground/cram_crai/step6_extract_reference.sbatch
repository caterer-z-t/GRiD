#!/bin/bash
#SBATCH --job-name=extract_ref_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-r REFERENCE_FA] [-b BED_FILE] [-o OUTPUT_DIR] [-f OUTPUT_PREFIX] [-h]"
    echo ""
    echo "  -r   Path to reference genome FASTA file (e.g., hs37d5.fa)"
    echo "  -b   Path to BED file with regions to extract (e.g., ref_hg37.bed)"
    echo "  -o   Output directory"
    echo "  -f   Output file prefix (default: ref_lpa)"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "r:b:o:f:h" opt; do
    case $opt in
        r) REFERENCE_FA="$OPTARG" ;;
        b) BED_FILE="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        f) OUTPUT_PREFIX="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Set default output prefix if not provided
OUTPUT_PREFIX="${OUTPUT_PREFIX:-ref_lpa}"

require_args REFERENCE_FA BED_FILE OUTPUT_DIR

export REFERENCE_FA BED_FILE OUTPUT_DIR OUTPUT_PREFIX

# -------------------------------
# Modules
# -------------------------------
module purge
module load bedtools/2.30.0

# -------------------------------
# Paths & Settings
# -------------------------------
require_file "$REFERENCE_FA"
require_file "$BED_FILE"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

OUTPUT_FASTA="${OUTPUT_DIR}/${OUTPUT_PREFIX}.fasta"

# -------------------------------
# Run Analysis
# -------------------------------
echo "Starting reference extraction..." > $OUTPUT_DIR/extract_reference.log 2>&1
echo "Reference: $REFERENCE_FA" >> $OUTPUT_DIR/extract_reference.log 2>&1
echo "BED file: $BED_FILE" >> $OUTPUT_DIR/extract_reference.log 2>&1
echo "Output: $OUTPUT_FASTA" >> $OUTPUT_DIR/extract_reference.log 2>&1
echo "" >> $OUTPUT_DIR/extract_reference.log 2>&1

# Extract reference sequence from specified regions
bedtools getfasta \
    -fi $REFERENCE_FA \
    -bed $BED_FILE \
    > $OUTPUT_FASTA 2>> $OUTPUT_DIR/extract_reference.log

if [ $? -eq 0 ]; then
    echo "Reference extraction completed successfully." >> $OUTPUT_DIR/extract_reference.log 2>&1
    echo "Output file: $OUTPUT_FASTA" >> $OUTPUT_DIR/extract_reference.log 2>&1
    echo "" >> $OUTPUT_DIR/extract_reference.log 2>&1
    echo "First few lines of output:" >> $OUTPUT_DIR/extract_reference.log 2>&1
    head -n 10 $OUTPUT_FASTA >> $OUTPUT_DIR/extract_reference.log 2>&1
else
    echo "ERROR: Reference extraction failed!" >> $OUTPUT_DIR/extract_reference.log 2>&1
    exit 1
fi