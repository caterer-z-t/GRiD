#!/bin/bash
#SBATCH --job-name=SOL_PIPELINE
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# ============================
#        USAGE
# sbatch sol_pipeline.sbatch        # runs steps 1–5 (default)
# sbatch sol_pipeline.sbatch 5      # runs only step 5
# sbatch sol_pipeline.sbatch 3-5    # runs steps 3 through 5
# ============================

# Parse input argument
RANGE=${1:-1-9}  # Default: full pipeline if nothing provided

# Convert argument to start and end numbers
if [[ "$RANGE" =~ ^([0-9]+)-([0-9]+)$ ]]; then
    START_STEP=${BASH_REMATCH[1]}
    END_STEP=${BASH_REMATCH[2]}
elif [[ "$RANGE" =~ ^[0-9]+$ ]]; then
    START_STEP=$RANGE
    END_STEP=$RANGE
elif [[ -z "$RANGE" ]]; then
    START_STEP=1
    END_STEP=5
else
    echo "❌ Invalid input: '$RANGE'"
    echo "Usage examples:"
    echo "  sbatch sol_pipeline.sbatch        # run all steps"
    echo "  sbatch sol_pipeline.sbatch 5      # run only step 5"
    echo "  sbatch sol_pipeline.sbatch 3-5    # run steps 3–5"
    exit 1
fi

# Helper function to determine if a step should run
run_step() {
    local step_num=$1
    if (( step_num >= START_STEP && step_num <= END_STEP )); then
        return 0  # Run
    else
        return 1  # Skip
    fi
}

echo "-------------------------------------------------"
echo " Running SOL pipeline steps $START_STEP through $END_STEP"
echo "-------------------------------------------------"

# ------------------------------
# Environment Setup
# ------------------------------
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate lpa_vntr

export WORK=/work/users/c/a/catererz
export USERS=/users/c/a/catererz
export HOME=/nas/longleaf/home/catererz
export EPI=/proj/epi
export SOL=$EPI/Genetic_Data_Center/SOL
export CRAM=$SOL/cram
export LPA=$HOME/epi/lpa
export GRID=$HOME/epi/GRiD
export SOFTWARE=$HOME/epi/software

# Name for working directory for this run
WORK_DIR=$WORK/LPA_VNTR_PIPELINE_SOFTWARE

mkdir -p $WORK_DIR
cd $WORK_DIR

REGIONS_FILE="$GRID/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

# ------------------------------
# Step 1: Ensure CRAM Indexes
# ------------------------------
if run_step 1; then
    echo "=== Step 1: Ensuring CRAM Indexes ==="
    for cram_file in $CRAM/*.cram; do
        if [ ! -f "${cram_file}.crai" ]; then
            echo "Index file for $cram_file not found. Creating index..." >> ensure_crai.log 2>&1
            grid ensure-crai --cram-file "$cram_file" >> ensure_crai.log 2>&1
        else
            echo "Index file for $cram_file already exists. Skipping..." >> ensure_crai.log 2>&1
        fi
    done
fi

# ------------------------------
# Step 2: Count Reads
# ------------------------------
if run_step 2; then
    echo "=== Step 2: Counting Reads ==="
    grid count-reads \
        --cram-dir $CRAM \
        --output-file count-reads.tsv \
        --ref-fasta $GRID/files/hg38.fa \
        --chrom chr$CHR \
        --start $START \
        --end $END \
        --config $GRID/grid/config.yaml \
        --threads 4 \
        >> count_reads.log 2>&1
fi

# ------------------------------
# Step 3: Run Mosdepth
# ------------------------------
if run_step 3; then
    echo "=== Step 3: Running Mosdepth ==="
    module load mosdepth
    grid mosdepth \
        --cram-dir $CRAM \
        --output-file mosdepth.tsv \
        --ref-fasta $GRID/files/hg38.fa \
        --chrom chr$CHR \
        --start $START \
        --end $END \
        --region-name LPA \
        --work-dir $WORK_DIR/mosdepth \
        --by 1000 \
        --fast \
        --threads 4 \
        >> mosdepth.log 2>&1
fi

# ------------------------------
# Step 4: Normalize Mosdepth
# ------------------------------
if run_step 4; then
    echo "=== Step 4: Normalizing Mosdepth ==="
    grid normalize-mosdepth \
        --mosdepth-dir $WORK_DIR/mosdepth \
        --repeat-mask $GRID/files/repeat_mask_list.hg38.ucsc_bed \
        --chrom chr$CHR \
        --start $START \
        --end $END \
        --threads 4 \
        --min-depth 20 \
        --max-depth 100 \
        --top-frac 0.1 \
        --output-file normalized_mosdepth.txt.gz \
        >> normalize_mosdepth.log 2>&1
fi

# ------------------------------
# Step 5: Find Neighbors
# ------------------------------
if run_step 5; then
    echo "=== Step 5: Finding Neighbors ==="
    grid find-neighbors \
        --input-file normalized_mosdepth.txt.gz \
        --output-file neighbors.txt \
        --zmax 2.0 \
        --n-neighbors 500 \
        --sigma2-max 1000.0 \
        >> find_neighbors.log 2>&1
fi

# ------------------------------
# Step 6: Extract LPA Reference
# ------------------------------
if run_step 6; then
    echo "=== Step 6: Extracting LPA Reference ==="
    grid extract-reference \
        --reference-fa $GRID/files/hg38.fa \
        --bed-file $GRID/files/ref_hg38.bed \
        --output-dir $WORK_DIR \
        --output-prefix lpa_reference \
        >> extract_lpa_reference.log 2>&1
fi

# ------------------------------
# Step 7: Realign Reads to LPA
# ------------------------------
if run_step 7; then
    echo "=== Step 7: Realigning Reads to LPA ==="
    grid lpa-realign \
        --cram-dir $CRAM \
        --output-file lpa_realigned.tsv \
        --ref-fasta $GRID/files/hg38.fa \
        --lpa-ref-fasta lpa_reference.fasta \
        --positions-file $GRID/files/hardcoded_positions.txt \
        --genome-build hg38 \
        --chrom chr$CHR \
        --start $START \
        --end $END \
        --threads 4 \
        >> lpa_realign.log 2>&1
fi


# ------------------------------
# Step 8: Estimate LPA DIPCN
# ------------------------------
if run_step 8; then
    echo "=== Step 8: Estimating LPA DIP CN ==="
    grid compute-dipcn \
        --count-file lpa_realigned.tsv \
        --neighbor-file neighbors.zMax2.0.txt.gz \
        --output-prefix lpa_dipcn \
        --n-neighbors 500 \
        >> compute_dipcn.log 2>&1
fi

# ------------------------------
# Step 9: Compute LPA KIV-2 CN
# ------------------------------
if run_step 9; then
    echo "=== Step 9: Computing LPA KIV-2 CN ==="
    grid estimate-kiv \
        --exon1a lpa_dipcn.exon1A.dipCN.txt \
        --exon1b lpa_dipcn.exon1B.dipCN.txt \
        --output lpa_kiv2_cn.txt \
        --format txt \
        >> compute_kiv2_cn.log 2>&1
fi


echo "Pipeline complete for steps $START_STEP–$END_STEP."
