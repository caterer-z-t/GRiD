#!/bin/bash
#SBATCH --job-name=realign_lpa
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-c CRAM_DIR] [-r REFERENCE_FA] [-f REF_FASTA] [-o OUTPUT_FILE] [-C CHR] [-s START] [-e END] [-h]"
    echo ""
    echo "  -p   Python realignment script"
    echo "  -c   Directory containing CRAM/BAM files"
    echo "  -r   Reference genome for CRAM files (e.g., hg38.fa)"
    echo "  -f   LPA reference FASTA (output from extract_reference)"
    echo "  -o   Output file path"
    echo "  -C   Chromosome (default: 6)"
    echo "  -s   Start position (default: 161032032)"
    echo "  -e   End position (default: 161068032)"
    echo "  -h   Show this help message"
    exit 1
}

# Default values for LPA region
CHR="6"
START="161032032"
END="161068032"

# Parse command line arguments
while getopts "p:c:r:f:o:C:s:e:h" opt; do
    case $opt in
        p) REALIGN_PY="$OPTARG" ;;
        c) CRAM_DIR="$OPTARG" ;;
        r) REFERENCE_FA="$OPTARG" ;;
        f) LPA_REF_FASTA="$OPTARG" ;;
        o) OUTPUT_FILE="$OPTARG" ;;
        C) CHR="$OPTARG" ;;
        s) START="$OPTARG" ;;
        e) END="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args REALIGN_PY CRAM_DIR REFERENCE_FA LPA_REF_FASTA OUTPUT_FILE

export REALIGN_PY CRAM_DIR REFERENCE_FA LPA_REF_FASTA OUTPUT_FILE CHR START END

# -------------------------------
# Modules
# -------------------------------
module purge
module load samtools/1.17

# Load Anaconda and activate environment
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate lpa_vntr

# -------------------------------
# Paths & Settings
# -------------------------------
require_file "$REALIGN_PY"
require_file "$REFERENCE_FA"
require_file "$LPA_REF_FASTA"
require_dir "$CRAM_DIR"

# Create output directory if needed
OUTPUT_DIR=$(dirname "$OUTPUT_FILE")
mkdir -p "$OUTPUT_DIR"

# Region to extract
REGION="${CHR}:${START}-${END}"

# Log file
LOG_FILE="${OUTPUT_DIR}/realign.log"

# -------------------------------
# Run Analysis
# -------------------------------
echo "Starting LPA realignment analysis..." > $LOG_FILE
echo "CRAM directory: $CRAM_DIR" >> $LOG_FILE
echo "Reference genome: $REFERENCE_FA" >> $LOG_FILE
echo "LPA reference: $LPA_REF_FASTA" >> $LOG_FILE
echo "Region: $REGION" >> $LOG_FILE
echo "Output: $OUTPUT_FILE" >> $LOG_FILE
echo "" >> $LOG_FILE

start_time=$(date +%s)

# Remove output file if it exists
rm -f "$OUTPUT_FILE"

# Find all CRAM/BAM files
CRAM_FILES=()
for ext in cram bam; do
    while IFS= read -r -d '' file; do
        CRAM_FILES+=("$file")
    done < <(find "$CRAM_DIR" -name "*.${ext}" -print0)
done

echo "Found ${#CRAM_FILES[@]} CRAM/BAM files" >> $LOG_FILE
echo "" >> $LOG_FILE

# Process each file
processed=0
failed=0

for CRAM_FILE in "${CRAM_FILES[@]}"; do
    # Extract sample ID from filename
    SAMPLE_ID=$(basename "$CRAM_FILE" | sed -E 's/\.(cram|bam)$//')
    
    # Check if index exists
    if [[ -f "${CRAM_FILE}.crai" ]] || [[ -f "${CRAM_FILE}.bai" ]]; then
        echo "Processing: $SAMPLE_ID" >> $LOG_FILE
        
        # Run samtools view piped to realignment script
        # Flags: 99=proper pair R1, 147=proper pair R2 reverse, 83=proper pair R1 reverse, 163=proper pair R2
        samtools view \
            --reference "$REFERENCE_FA" \
            "$CRAM_FILE" \
            "$REGION" \
            --verbosity 0 \
            2>> $LOG_FILE \
            | awk '($2==99||$2==147||$2==83||$2==163) {print $1,$4,$6,$10,$11}' \
            | sort -k1,1 -k2,2n \
            | python "$REALIGN_PY" \
                --reference "$LPA_REF_FASTA" \
                --sample_id "$SAMPLE_ID" \
                >> "$OUTPUT_FILE" 2>> $LOG_FILE
        
        if [ $? -eq 0 ]; then
            ((processed++))
        else
            echo "WARNING: Failed to process $SAMPLE_ID" >> $LOG_FILE
            ((failed++))
        fi
    else
        echo "WARNING: No index found for $CRAM_FILE, skipping" >> $LOG_FILE
        ((failed++))
    fi
done

end_time=$(date +%s)
elapsed=$((end_time - start_time))

echo "" >> $LOG_FILE
echo "Realignment completed!" >> $LOG_FILE
echo "Successfully processed: $processed files" >> $LOG_FILE
echo "Failed/skipped: $failed files" >> $LOG_FILE
echo "Time elapsed: ${elapsed} seconds" >> $LOG_FILE
echo "" >> $LOG_FILE
echo "First few lines of output:" >> $LOG_FILE
head -n 20 "$OUTPUT_FILE" >> $LOG_FILE
