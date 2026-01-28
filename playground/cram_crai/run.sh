#!/bin/bash

# Step 1, move files from google bucket to local (aka longleaf)
# sbatch step1_gcloud_copy_files.sbatch -d gs://fc-secure-708ee835-1980-4752-941f-f42b5ce2f58b/submissions/

# Step 2, count reads in VNTR region for each cram file
REGIONS_FILE="$GRID/files/734_possible_coding_vntr_regions.IBD2R_gt_0.25.uniq.txt"
LPA_REGION=$(awk '$7=="LPA" {print $1,$2,$3}' "$REGIONS_FILE")
CHR=$(echo "$LPA_REGION" | awk '{print $1}')
START=$(echo "$LPA_REGION" | awk '{print $2}')
END=$(echo "$LPA_REGION" | awk '{print $3}')

sbatch step2_count_read.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_read_counts/read_count.txt \
   -r $GRID/files/hg38.fa \
   -C $CHR \
   -s $START \
   -e $END \
   -w $WORK/lpa_read_counts

# Step 3, run mosdepth to get coverage across VNTR region for each cram file
sbatch step3_mosdepth.sbatch \
   -c $SOL/cram \
   -o $WORK/lpa_mosdepth/mosdepth_summary.txt \
   -r $GRID/files/hg38.fa \
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

# Step 6, Extract reference sequences for VNTR region
# Step LPA-1, extract reference sequence for LPA region
sbatch step6_extract_reference.sbatch \
   -r $GRID/files/hg38.fa \
   -b $GRID/files/ref_hg38.bed \
   -o $WORK/lpa_extract_reference \
   -f ref_lpa_hg38

# Step 7, realign reads to extracted reference sequence
# Step LPA-2, realign reads to LPA reference sequence
sbatch step7_realign.sbatch \
   -p $GRID/playground/cram_crai/src/realign_lpa.py \
   -c $SOL/cram \
   -r $GRID/files/hg38.fa \
   -f $WORK/lpa_extract_reference/ref_lpa_hg38.fasta \
   -o $WORK/lpa_realign/realign_results.txt \
   -P $GRID/files/hardcoded_positions.txt \
   -g hg38 \
   -C $CHR \
   -s $START \
   -e $END


# Step 8, compute diploid CN estimates
# Step LPA-3, compute diploid CN estimates for LPA region
sbatch step8_compute_dipCN.sbatch \
   -p $GRID/playground/cram_crai/src/compute_dipCN.py \
   -c $WORK/lpa_realign/realign_results.txt \
   -n $WORK/lpa_find_neighbors/LPA_neighbors.zMax1.0.txt.gz \
   -o $WORK/lpa_compute_dipCN/LPA \
   -N 200

# Step 9, estimate KIV CN
# Step LPA-4, estimate KIV CN for LPA region
sbatch step9_estimate_KIV_CN.sbatch \
   -p $GRID/playground/cram_crai/src/estimate_KIV_CN.py \
   -a $WORK/lpa_compute_dipCN/LPA.exon1A.dipCN.txt \
   -b $WORK/lpa_compute_dipCN/LPA.exon1B.dipCN.txt \
   -o $WORK/lpa_estimate_kiv_cn/LPA_KIV2_copy_numbers.tsv \
   -f tsv

# Step 10, plot KIV2 copy number vs Lp(a) association
sbatch step10_plot_KIVCN_and_lpa.sbatch \
   -p $GRID/playground/cram_crai/src/plot_lpa_ass.py \
   -k $WORK/lpa_estimate_kiv_cn/LPA_KIV2_copy_numbers.tsv \
   -l $GRID/files/diploid_calls.txt \
   -o $WORK/lpa_plots/LPA_association