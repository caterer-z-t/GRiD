#!/bin/bash
#SBATCH --job-name=normalize_neighbors
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# --- Activate environment ---
module purge
module load anaconda
source $(conda info --base)/etc/profile.d/conda.sh 
eval "$(conda shell.bash hook)" 
conda activate lpa_vntr

# --- File paths (hg38 version) ---
COUNT_READ_OUTPUT_PATH="$WORK/lpa_read_counts/VNTR_counts.txt"
FIND_NEIGHBOR_OUTPUT_PATH="$WORK/lpa_neighbors/LPA_find_neighbor.zMax2.0.txt.gz"
OUTPUT_PATH="$WORK/lpa_normalize_neighbors"
IBD2R_path="$LPA/playground/cram_crai/step4_normalize_mosdepth/repeat_mask_list.hg38.ucsc_bed"
PREFIX="hg38_20250919"

# make sure output directory exists
mkdir -p $OUTPUT_PATH
cd $OUTPUT_PATH || { echo "ERROR: Cannot change to output directory: ${OUTPUT_PATH}"; exit 1; }

COL=1

# Loop through VNTR regions in IBD2R file
while read LINE; do
    COL=$((COL+1))
    ARR=($LINE)

    CHR=${ARR[0]}
    BP_START=${ARR[1]}
    BP_END=${ARR[2]}
    VNTR=${PREFIX}_${CHR}_${BP_START}_${BP_END}

    echo "Processing column $COL -> $VNTR"

    N_NBR=300

    awk -v Nnbr=$N_NBR -v col=$COL -v chr=$CHR '
    ARGIND==1 {reads[$1]=$col}
    ARGIND==2 {
      sum=0; num=0;
      for (i=1; i<=Nnbr; i++) {
          sum+=reads[$(3*i)]/$(3*i+1);
          num++;
      }
      norm_reads = reads[$1]/$2/(sum/num);
      print $1, norm_reads
    }' \
    $COUNT_READ_OUTPUT_PATH \
    <(zcat $FIND_NEIGHBOR_OUTPUT_PATH) \
    > $OUTPUT_PATH/$VNTR.dipCN.txt

done < <(tail -n +2 $IBD2R_path)
