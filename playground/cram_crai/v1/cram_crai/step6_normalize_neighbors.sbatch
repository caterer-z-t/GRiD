#!/bin/bash
#SBATCH --job-name=normalize_neighbors
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=general
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

source "$GRID/playground/cram_crai/common.sh"

usage(){
    echo "Usage: $0 -c COUNT_FILE -n NEIGHBOR_FILE -o OUTPUT_DIR -r REGION_FILE -p PREFIX -C CHR -s START -e END [-N N_NEIGHBORS]"
    echo ""
    echo "  -c   VNTR read count file (from step2)"
    echo "  -n   Neighbor results file (from step5, .txt.gz)"
    echo "  -o   Output directory"
    echo "  -r   VNTR regions file (IBD2R bed-like file)"
    echo "  -p   Prefix for output filenames"
    echo "  -C   Chromosome of VNTR"
    echo "  -s   Start position of VNTR"
    echo "  -e   End position of VNTR"
    echo "  -N   Number of neighbors to average over (default: 300)"
    echo "  -h   Show this help message"
    exit 1
}

# -------------------------------
# Parse command line arguments
# -------------------------------
# Command line args
while getopts c:n:o:r:p:C:s:e:N: flag
do
    case "${flag}" in
        c) COUNT_READ_OUTPUT_PATH=${OPTARG};;
        n) FIND_NEIGHBOR_OUTPUT_PATH=${OPTARG};;
        o) OUTPUT_PATH=${OPTARG};;
        r) IBD2R_path=${OPTARG};;
        p) PREFIX=${OPTARG};;
        C) CHR=${OPTARG};;
        s) START=${OPTARG};;
        e) END=${OPTARG};;
        N) N_NBR=${OPTARG};;
    esac
done

require_args COUNT_READ_OUTPUT_PATH FIND_NEIGHBOR_OUTPUT_PATH OUTPUT_PATH IBD2R_path PREFIX CHR START END

export COUNT_READ_OUTPUT_PATH FIND_NEIGHBOR_OUTPUT_PATH OUTPUT_PATH IBD2R_path PREFIX CHR START END

# -------------------------------
# Modules
# -------------------------------
module purge

# -------------------------------
# Paths & Settings
# -------------------------------
require_file "$COUNT_READ_OUTPUT_PATH"
require_file "$FIND_NEIGHBOR_OUTPUT_PATH"
require_file "$IBD2R_path"

mkdir -p $OUTPUT_PATH
cd $OUTPUT_PATH || exit 1

COL=1

# lets print the head of each of the files
awk 'NR==1 {print; next} {sub(/\..*/,"",$1); print}' \
  $COUNT_READ_OUTPUT_PATH \
  > read_count.cleaned.txt

# echo "First 3 lines of read_count.cleaned.txt:"
# head -n 3 read_count.cleaned.txt 

# echo "FIND_NEIGHBOR_OUTPUT_PATH: $FIND_NEIGHBOR_OUTPUT_PATH"
# zcat $FIND_NEIGHBOR_OUTPUT_PATH | head -n 3

# echo "IBD2R_path: $IBD2R_path"
# awk -v chr="$CHR" -v start="$START" -v end="$END" '$1==chr && $2==start && $3==end' "$IBD2R_path" | head -n 3

# -------------------------------
# Run normalization
# -------------------------------
COL=1
N_NBR=${N_NBR:-300}   # default if not set

# Loop over matching VNTRs in the region file
awk -v chr="$CHR" -v start="$START" -v end="$END" '$1==chr && $2==start && $3==end' "$IBD2R_path" | \
while read CHR BP_START BP_END REST; do
    COL=$((COL+1))
    VNTR=${PREFIX}_${CHR}_${BP_START}_${BP_END}

    echo "Processing $VNTR using column $COL and $N_NBR neighbors..."

    awk -v Nnbr=$N_NBR -v col=$COL -v chr=$CHR '
    ARGIND==1 {
        sid=$1
        reads[sid]=$col
        next
    }
    ARGIND==2 {
        # Only process if this sample exists in read_count file
        if (!($1 in reads)) {
            next
        }
        
        # Only process if we have a valid read count
        if (reads[$1] == "" || reads[$1] == 0) {
            next
        }
        
        # Calculate neighbor average
        sum=0; num=0;
        for (i=1; i<=Nnbr; i++) {
            nbr_sample = $(3*i)
            nbr_dist = $(3*i+1)
            
            # Check if neighbor exists and has valid data
            if (nbr_sample in reads && reads[nbr_sample] != "" && reads[nbr_sample] != 0 && nbr_dist > 0) {
                sum += reads[nbr_sample] / nbr_dist
                num++
            }
        }
        
        # Only output if we have valid neighbors
        if (num > 0 && $2 > 0) {
            norm_reads = reads[$1] / $2 / (sum / num)
            print $1, norm_reads
        }
    }' \
    read_count.cleaned.txt \
    <(zcat "$FIND_NEIGHBOR_OUTPUT_PATH") \
    > "$OUTPUT_PATH/$VNTR.dipCN.txt"

done