#!/bin/bash

# LPA VNTR Pipeline - Dependency-aware submission
# This script submits all jobs with proper dependencies so they run sequentially
# Each job waits for the previous one to complete successfully before starting

# set -e  # Exit on error

echo "=========================================="
echo "LPA VNTR Pipeline - Dependency Submission"
echo "=========================================="
echo ""

WORK="/nas/longleaf/home/catererz/work/LPA_VNTR_PIPELINE_CODE"
mkdir -p $WORK

# Get region information
REGIONS_FILE="$GRID/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

echo "LPA Region: Chr${CHR}:${START}-${END}"
echo ""

# Step 1: Copy files from Google bucket
# echo "Submitting Step 1: Copy files from Google bucket..."
# JOB1=$(sbatch --parsable step1_gcloud_copy_files.sbatch -d gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/)
# echo "  Job ID: $JOB1"
# echo ""

# Step 2: Count reads in VNTR region
echo "Submitting Step 2: Count reads (depends on Step 1)..."
JOB2=$(sbatch --parsable step2_count_read.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_read_counts/read_count.txt \
   -r $GRID/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_read_counts)
echo "  Job ID: $JOB2"
echo ""

# Step 3: Run mosdepth
echo "Submitting Step 3: Mosdepth coverage (depends on Step 1)..."
JOB3=$(sbatch --parsable step3_mosdepth.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_mosdepth/mosdepth_summary.txt \
   -r $GRID/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_mosdepth)
echo "  Job ID: $JOB3"
echo ""

# Step 4: Normalize mosdepth
echo "Submitting Step 4: Normalize mosdepth (depends on Step 3)..."
JOB4=$(sbatch --parsable --dependency=afterok:$JOB3 step4_normalize_mosdepth.sbatch \
   -p $GRID/playground/cram_crai/src/normalize_mosdepth.py \
   -m $WORK/lpa_mosdepth \
   -d $WORK/lpa_norm_mosdepth \
   -r $GRID/files/repeat_mask_list.hg38.ucsc_bed \
   -c $CHR \
   -s $START \
   -e $END \
   -n 20 \
   -x 100)
echo "  Job ID: $JOB4"
echo ""

# Step 5: Find neighbors
echo "Submitting Step 5: Find neighbors (depends on Step 4)..."
JOB5=$(sbatch --parsable --dependency=afterok:$JOB4 step5_find_neighbors.sbatch \
   -p $GRID/playground/cram_crai/src/find_neighbors.py \
   -o $WORK/lpa_find_neighbors \
   -i $WORK/lpa_norm_mosdepth/normalized_output.txt.gz \
   -f LPA_neighbors \
   -z 1.0 \
   -n 5)
echo "  Job ID: $JOB5"
echo ""

# Step LPA-1: Extract reference
echo "Submitting Step LPA-1: Extract reference (no dependencies)..."
JOB_LPA1=$(sbatch --parsable step6_extract_reference.sbatch \
   -r $GRID/cram/hg38.fa \
   -b $GRID/files/ref_hg37.bed \
   -o $WORK/lpa_extract_reference \
   -f ref_lpa_hg38)
echo "  Job ID: $JOB_LPA1"
echo ""

# Step LPA-2: Realign
# Depends on Step 1 (CRAM files) AND Step LPA-1 (reference)
echo "Submitting Step LPA-2: Realign reads (depends on Step 1 and LPA-1)..."
JOB_LPA2=$(sbatch --parsable --dependency=afterok:$JOB_LPA1 step7_realign.sbatch \
   -p $GRID/playground/cram_crai/src/realign_lpa.py \
   -c $SOL/cram \
   -r $GRID/cram/hg38.fa \
   -f $WORK/lpa_extract_reference/ref_lpa_hg38.fasta \
   -o $WORK/lpa_realign/realign_results.txt \
   -C $CHR \
   -s $START \
   -e $END)
echo "  Job ID: $JOB_LPA2"
echo ""

# Step LPA-3: Compute diploid CN
# Depends on Step 5 (neighbors) AND Step LPA-2 (realignment)
echo "Submitting Step LPA-3: Compute diploid CN (depends on Step 5 and LPA-2)..."
JOB_LPA3=$(sbatch --parsable --dependency=afterok:$JOB5:$JOB_LPA2 step8_compute_dipCN.sbatch \
   -p $GRID/playground/cram_crai/src/compute_dipCN.py \
   -c $WORK/lpa_realign/realign_results.txt \
   -n $WORK/lpa_find_neighbors/LPA_neighbors.zMax1.0.txt.gz \
   -o $WORK/lpa_compute_dipCN/LPA \
   -N 200)
echo "  Job ID: $JOB_LPA3"
echo ""

# Step LPA-5: Estimate KIV CN
echo "Submitting Step LPA-5: Estimate KIV CN (depends on Step LPA-3)..."
JOB_LPA5=$(sbatch --parsable --dependency=afterok:$JOB_LPA3 step9_estimate_KIV_CN.sbatch \
   -p $GRID/playground/cram_crai/src/estimate_KIV_CN.py \
   -a $WORK/lpa_compute_dipCN/LPA.exon1A.dipCN.txt \
   -b $WORK/lpa_compute_dipCN/LPA.exon1B.dipCN.txt \
   -o $WORK/lpa_estimate_kiv_cn/LPA_KIV2_copy_numbers.tsv \
   -f tsv)
echo "  Job ID: $JOB_LPA5"
echo ""

echo "=========================================="
echo "All jobs submitted successfully!"
echo "=========================================="
echo ""
echo "Pipeline dependency chain:"
echo "  Step 1 (Copy files)        → $JOB1"
echo "  ├─ Step 2 (Count reads)    → $JOB2"
echo "  ├─ Step 3 (Mosdepth)       → $JOB3"
echo "  │  └─ Step 4 (Normalize)   → $JOB4"
echo "  │     └─ Step 5 (Neighbors)→ $JOB5"
echo "  │        └─┐"
echo "  │          │"
echo "  └─────────┐│"
echo "            ││"
echo "  LPA-1 (Extract ref)        → $JOB_LPA1"
echo "  └─┬───────┘│"
echo "       │     │"
echo "  LPA-2 (Realign)            → $JOB_LPA2"
echo "       └─────┤"
echo "             │"
echo "  LPA-3 (Diploid CN)         → $JOB_LPA3"
echo "   └─┬───────┘"
echo "     │"
echo "  LPA-5 (Estimate KIV CN)    → $JOB_LPA5 (FINAL)"
echo ""
echo "Monitor job status with: squeue -u $USER"
echo "View job details with: scontrol show job <JOB_ID>"
echo "Cancel all jobs with: scancel $JOB1 $JOB2 $JOB3 $JOB4 $JOB5 $JOB_LPA1 $JOB_LPA2 $JOB_LPA3 $JOB_LPA5"
echo ""