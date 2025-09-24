#!/bin/bash

# Step 1, move files from google bucket to local (aka longleaf)
sbatch step1_gcloud_copy_files.sbatch -d gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/

# Step 2, count reads in VNTR region for each cram file
REGIONS_FILE="$LPA/playground/cram_crai/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

sbatch step2_count_read.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_read_counts/read_count.txt \
   -r $LPA/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_read_counts

# Step 3, run mosdepth to get coverage across VNTR region for each cram file
sbatch step3_mosdepth.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_mosdepth/mosdepth_summary.txt \
   -r $LPA/cram/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_mosdepth

# Step 4, normalize mosdepth coverage values
sbatch step4_normalize_mosdepth.sbatch \
   -p $LPA/playground/cram_crai/src/normalize_mosdepth.py \
   -m $WORK/lpa_mosdepth \
   -d $WORK/lpa_norm_mosdepth \
   -r $LPA/playground/cram_crai/files/repeat_mask_list.hg38.ucsc_bed \
   -c $CHR \
   -s $START \
   -e $END \
   -n 20 \
   -x 100 

# Step 5, Find Neighbors
sbatch step5_find_neighbors.sbatch \
   -p $LPA/playground/cram_crai/src/find_neighbors.py \
   -o $WORK/lpa_find_neighbors \
   -i $WORK/lpa_norm_mosdepth_multi/normalized_output.txt.gz \
   -f LPA_neighbors \
   -z 1.0 \
   -n 5