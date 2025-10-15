#!/bin/bash
#SBATCH --job-name=normalize_mosdepth
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=10:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage() {
    echo "Usage: $0 [-p PYTHON_SCRIPT] [-m MOSDEPTH_DIR] [-d OUTPUT_DIR] [-r REPEAT_MASK] [-c CHROMOSOME] [-s START] [-e END] [-n MIN_DEPTH] [-x MAX_DEPTH] [-h]"
    echo ""
    echo "  -p   Python script to run"
    echo "  -m   Directory containing mosdepth files"
    echo "  -d   Output directory"
    echo "  -c   Chromosome to analyze"
    echo "  -s   Start position"
    echo "  -e   End position"
    echo "  -n   Minimum depth threshold"
    echo "  -x   Maximum depth threshold"
    echo "  -r   Repeat mask BED file"
    echo "  -h   Show this help message"
    exit 1
}

# Parse command line arguments
while getopts "p:m:d:r:c:s:e:n:x:h" opt; do
    case $opt in
        p) PYTHON_SCRIPT="$OPTARG" ;;
        m) MOSDEPTH_DIR="$OPTARG" ;;
        d) OUT_DIR="$OPTARG" ;;
        r) REPEAT_MASK="$OPTARG" ;;
        c) CHROMOSOME="$OPTARG" ;;
        s) START="$OPTARG" ;;
        e) END="$OPTARG" ;;
        n) MIN_DEPTH="$OPTARG" ;;
        x) MAX_DEPTH="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

require_args PYTHON_SCRIPT MOSDEPTH_DIR OUT_DIR REPEAT_MASK CHROMOSOME START END MIN_DEPTH MAX_DEPTH

export PYTHON_SCRIPT MOSDEPTH_DIR OUT_DIR REPEAT_MASK CHROMOSOME START END MIN_DEPTH MAX_DEPTH

# -------------------------------
# Modules
# -------------------------------
module purge

# Load Anaconda and activate environment
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate lpa_vntr

# -------------------------------
# Paths & Settings
# -------------------------------
require_file "$PYTHON_SCRIPT"
require_file "$REPEAT_MASK"

require_dir "$MOSDEPTH_DIR"

mkdir -p $OUT_DIR

# -------------------------------
# Run Python Normalization
# -------------------------------

# Run the normalization script with correct arguments
# python $PYTHON_SCRIPT \
#     --mosdepth_dir $MOSDEPTH_DIR \
#     --repeat_mask $REPEAT_MASK \
#     --chromosome $CHROMOSOME \
#     --start $START \
#     --end $END \
#     --min_depth $MIN_DEPTH \
#     --max_depth $MAX_DEPTH \
#     --output_file $OUT_DIR/normalized_output.txt.gz \
#     > $OUT_DIR/normalize_mosdepth_python.log 2>&1

# zcat $OUT_DIR/normalized_output.txt.gz | head -20 | column -t

# echo "Normalization complete. Output written to $OUT_DIR/normalized_output.txt.gz"

# Run multi_threaded version
MULTI_THREADED_SCRIPT="$(dirname $PYTHON_SCRIPT)/multi_threaded_norm_mosdepth.py"

require_file "$MULTI_THREADED_SCRIPT"

python $MULTI_THREADED_SCRIPT \
    --mosdepth_dir $MOSDEPTH_DIR \
    --repeat_mask $REPEAT_MASK \
    --chromosome $CHROMOSOME \
    --start $START \
    --end $END \
    --threads ${SLURM_CPUS_PER_TASK:-1} \
    --min_depth $MIN_DEPTH \
    --max_depth $MAX_DEPTH \
    --output_file $OUT_DIR/normalized_output.txt.gz \
    > $OUT_DIR/normalize_mosdepth_python.log 2>&1

zcat $OUT_DIR/normalized_output.txt.gz | head -20 | column -t

echo "Normalization complete. Output written to $OUT_DIR/normalized_output.txt.gz"