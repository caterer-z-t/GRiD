#!/bin/bash
#SBATCH --job-name=ilash_allchr
#SBATCH --output=slurm/ilash.out
#SBATCH --error=slurm/ilash.err
#SBATCH --time=11-00:00:00
#SBATCH --mem=300G
#SBATCH --cpus-per-task=32
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

# ===============================================================
# iLASH Pipeline: Merge → Normalize → PLINK → iLASH
# Author: Zach Caterer (revised with GPT-5)
# ===============================================================

set -euo pipefail
echo "Job started on $(hostname) at $(date)"
echo "Working directory: $(pwd)"

# ---------------------------
# Environment setup
# ---------------------------
module purge
module load anaconda
module load bcftools
module load samtools
module load plink/1.90b6.21
module load ilash

source $(conda info --base)/etc/profile.d/conda.sh
conda activate gnomix

# ---------------------------
# Variables
# ---------------------------
CHROM="6"
WORKDIR="$WORK/ibd"
GNOMIX_OUTDIR="$WORK/gnomix/output"
SOFTWAREDIR="$SOFTWARE/iLASH"

mkdir -p "$WORKDIR"
cd "$WORKDIR"

# ---------------------------
# Determine ILASH executable
# ---------------------------
if command -v ilash &> /dev/null; then
    ILASH="$(which ilash)"
    echo "Using iLASH from module: $ILASH"
else
    ILASH="$SOFTWAREDIR/ilash"
    echo "Using iLASH from software dir: $ILASH"
fi

# ---------------------------
# Step 1: Merge GnoFix VCF chunks
# ---------------------------
echo ""
echo "=== Step 1: Merging VCF chunks for chr${CHROM} ==="
MERGED_VCF="phased_chr${CHROM}.vcf.gz"

if [[ -f "$MERGED_VCF" && -f "${MERGED_VCF}.tbi" ]]; then
    echo "Merged VCF already exists."
else
    VCF_LIST="${WORKDIR}/vcf_list_chr${CHROM}.txt"
    find "${GNOMIX_OUTDIR}" -maxdepth 1 -type f -name "gnofix_output_chr${CHROM}_*.vcf" | sort > "$VCF_LIST"

    NUM_FILES=$(wc -l < "$VCF_LIST" || echo 0)
    if [[ "$NUM_FILES" -eq 0 ]]; then
        echo "ERROR: No GnoFix VCFs found for chr${CHROM}!"
        exit 1
    fi

    echo "Found $NUM_FILES VCF chunks."

    # Compress + index if needed
    while read -r vcf; do
        if [[ ! -f "${vcf}.gz" ]]; then
            echo "Compressing $(basename "$vcf")..."
            bgzip -c "$vcf" > "${vcf}.gz"
            tabix -p vcf "${vcf}.gz"
        fi
    done < "$VCF_LIST"

    # Merge by variant positions
    echo "Merging..."
    sed 's/$/.gz/' "$VCF_LIST" > "${VCF_LIST}.gz"
    bcftools concat --file-list "${VCF_LIST}.gz" -Oz -o "$MERGED_VCF"
    tabix -p vcf "$MERGED_VCF"
fi

echo "VCF merged successfully."

# ---------------------------
# Step 2: Normalize + clean
# ---------------------------
echo ""
echo "=== Step 2: Normalizing and cleaning ==="
NORM_VCF="phased_chr${CHROM}_norm.vcf.gz"
DEDUP_VCF="phased_chr${CHROM}_dedup.vcf.gz"
FINAL_VCF="phased_chr${CHROM}_final.vcf.gz"

if [[ ! -f "$FINAL_VCF" ]]; then
    echo "Splitting multi-allelics..."
    bcftools norm -m -any "$MERGED_VCF" -Oz -o "$NORM_VCF" --threads $SLURM_CPUS_PER_TASK
    tabix -p vcf "$NORM_VCF"

    echo "Removing duplicates..."
    bcftools norm -d any "$NORM_VCF" -Oz -o "$DEDUP_VCF"
    tabix -p vcf "$DEDUP_VCF"

    echo "Assigning variant IDs and fixing chr labels..."
    bcftools annotate --set-id +'%CHROM\_%POS\_%REF\_%ALT' "$DEDUP_VCF" -Oz -o tmp1.vcf.gz
    bcftools annotate --rename-chrs <(echo -e "chr${CHROM}\t${CHROM}") -Oz -o "$FINAL_VCF" tmp1.vcf.gz
    tabix -p vcf "$FINAL_VCF"
    rm tmp1.vcf.gz
else
    echo "Normalized VCF already exists."
fi

VARIANT_COUNT=$(bcftools index -n "$FINAL_VCF")
SAMPLE_COUNT=$(bcftools query -l "$FINAL_VCF" | wc -l)
echo "Final VCF: $VARIANT_COUNT variants, $SAMPLE_COUNT samples"

# ---------------------------
# Step 3: PLINK conversion
# ---------------------------
echo ""
echo "=== Step 3: Converting to PLINK ==="
BASENAME="phased_chr${CHROM}"

if [[ ! -f "${BASENAME}.ped" || ! -f "${BASENAME}.map" ]]; then
    plink --vcf "$FINAL_VCF" \
          --recode \
          --out "$BASENAME" \
          --allow-extra-chr \
          --keep-allele-order \
          --double-id \
          --threads $SLURM_CPUS_PER_TASK
fi

# ---------------------------
# Step 4: Genetic map scaling
# ---------------------------
echo ""
echo "=== Step 4: Updating genetic map ==="

if [[ ! -f "${BASENAME}.map.orig" ]]; then
    cp "${BASENAME}.map" "${BASENAME}.map.orig"
fi

CHR_LENGTH=$(bcftools query -f '%POS\n' "$FINAL_VCF" | tail -1)
CHR_CM=180
awk -v len="$CHR_LENGTH" -v cm="$CHR_CM" '{printf "%s\t%s\t%.6f\t%s\n", $1, $2, ($4/len)*cm, $4}' \
    "${BASENAME}.map.orig" > "${BASENAME}.map"

# ---------------------------
# Step 5: Validation
# ---------------------------
echo ""
echo "=== Step 5: Validating PLINK files ==="
MAP_COUNT=$(wc -l < "${BASENAME}.map")
EXPECTED_FIELDS=$((MAP_COUNT * 2 + 6))
FIRST_SAMPLE_FIELDS=$(head -1 "${BASENAME}.ped" | awk '{print NF}')

if [[ "$FIRST_SAMPLE_FIELDS" -ne "$EXPECTED_FIELDS" ]]; then
    echo "ERROR: PED field count mismatch!"
    exit 1
fi

echo "✓ PLINK structure valid."

# ---------------------------
# Step 6: iLASH config + run
# ---------------------------
echo ""
echo "=== Step 6: Running iLASH ==="

CONFIG="ilash_config_chr${CHROM}.txt"
cat > "$CONFIG" <<EOF
map ${BASENAME}.map
ped ${BASENAME}.ped
output output_chr${CHROM}.match
auto_slice 1
slice_length 2.9
cm_overlap 1.0
min_length 2.9
shingle_size 20
shingle_overlap 0
perm_count 20
bucket_count 5
interest_threshold 0.70
match_threshold 0.99
minhash_threshold 55
max_thread ${SLURM_CPUS_PER_TASK}
EOF

echo "Running iLASH for chr${CHROM}..."
set +e
$ILASH "$CONFIG" 2>&1 | tee "ilash_run_chr${CHROM}.log"
EXIT_CODE=$?
set -e

if [[ "$EXIT_CODE" -ne 0 ]]; then
    echo "❌ iLASH failed for chr${CHROM} (exit code $EXIT_CODE)"
    tail -50 "ilash_run_chr${CHROM}.log"
    exit 1
else
    echo "✅ iLASH completed successfully for chr${CHROM}"
fi

# ---------------------------
# Step 7: Summary
# ---------------------------
OUTPUT_FILE="output_chr${CHROM}.match"
if [[ -f "$OUTPUT_FILE" ]]; then
    MATCH_COUNT=$(wc -l < "$OUTPUT_FILE")
    echo "IBD segments found: $MATCH_COUNT"
    awk '{sum+=$10; count++} END {
        if (count>0) print "Average segment length: " sum/count " cM"
    }' "$OUTPUT_FILE"
else
    echo "WARNING: iLASH output missing!"
fi

echo "=== Chromosome ${CHROM} completed at $(date) ==="
