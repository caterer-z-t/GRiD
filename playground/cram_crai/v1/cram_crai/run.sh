#!/bin/bash

# Step 1, move files from google bucket to local (aka longleaf)
sbatch step1_gcloud_copy_files.sbatch -d gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/

# Step 2, count reads in VNTR region for each cram file
REGIONS_FILE="$GRID/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

sbatch step2_count_read.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_read_counts/read_count.txt \
   -r $GRID/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_read_counts

# Step 3, run mosdepth to get coverage across VNTR region for each cram file
sbatch step3_mosdepth.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_mosdepth/mosdepth_summary.txt \
   -r $GRID/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_mosdepth

# Step 4, normalize mosdepth coverage values
sbatch step4_normalize_mosdepth.sbatch \
   -p $GRID/playground/cram_crai/src/normalize_mosdepth.py \
   -m $WORK/lpa_mosdepth \
   -d $WORK/lpa_norm_mosdepth \
   -r $GRID/files/repeat_mask_list.hg38.ucsc_bed \
   -c $CHR \
   -s $START \
   -e $END \
   -n 20 \
   -x 100 

# Step 5, Find Neighbors
sbatch step5_find_neighbors.sbatch \
   -p $GRID/playground/cram_crai/src/find_neighbors.py \
   -o $WORK/lpa_find_neighbors \
   -i $WORK/lpa_norm_mosdepth/normalized_output.txt.gz \
   -f LPA_neighbors \
   -z 1.0 \
   -n 5

# Step 6, Normalize neighbors
sbatch step6_normalize_neighbors.sbatch \
   -c $WORK/lpa_read_counts/read_count.txt \
   -n $WORK/lpa_find_neighbors/LPA_neighbors.zMax1.0.txt.gz \
   -o $WORK/lpa_normalize_neighbors \
   -r $GRID/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt \
   -p test \
   -C $CHR \
   -s $START \
   -e $END \
   -N 5

sbatch step7_compare_lpa_ass.sbatch \
    -p $GRID/playground/cram_crai/src/plot_vntr_lpa_dipCN_ass.py \
    -f $GRID/files/diploid_calls.txt \
    -d $WORK/lpa_normalize_neighbors/test_6_160605062_160647661.dipCN.txt \
    -o lpa_vntr_lpa_dipCN_ass \
    -O $WORK/lpa_vntr_lpa_dipCN_ass
