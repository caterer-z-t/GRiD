#!/bin/bash
#SBATCH --job-name=ilash_ch6
#SBATCH --output=slurm/ilash_ch6.out
#SBATCH --error=slurm/ilash_ch6.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=10
#SBATCH --mem=500g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
module load gcc/12.2.0
module load ilash/1.0.2
module load anaconda
# source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

# Define paths
WORKDIR=$WORK/ilash_ch6
OUTPUT_DIR=$LPA/ilash/ilash_ch6

# Define Test_File paths
FILEPREFIX=lpa_gene_phased
FILEPATHS=$LPA/ilash/convert_to_ped/
MAP=$FILEPATHS/$FILEPREFIX.map
PED=$FILEPATHS/$FILEPREFIX.ped
CONFIG_FILE=config.txt

# Make directories only if they don't exist
[ -d "$WORKDIR" ] || mkdir -p "$WORKDIR"
[ -d "$OUTPUT_DIR" ] || mkdir -p "$OUTPUT_DIR"

# remove all files in output directory except for the .sbatch file
find "$OUTPUT_DIR" -type f ! -name "*.sbatch" -delete

cd "$WORKDIR" || exit 1  # exit if cd fails

# Copy and decompress input files only if they don't already exist in $WORKDIR
[ -f "$(basename "$MAP")" ] || cp "$MAP" .
[ -f "$(basename "$PED")" ] || cp "$PED" .


# Create config
echo "⚙️ Creating iLASH config"
cat > "$CONFIG_FILE" <<EOF
map $MAP
ped $PED
output output.match
slice_size 350
step_size 350
perm_count 50
shingle_size 15
shingle_overlap 0
bucket_count 10
max_thread 10
match_threshold 0.95
interest_threshold 0.70
min_length 2.9
auto_slice 1
slice_length 2.9
cm_overlap 1
minhash_threshold 55
EOF

# Run iLASH
echo "🚀 Running iLASH"
ilash "$CONFIG_FILE"

mv -f "$WORKDIR"/* "$OUTPUT_DIR"

# Cleanup
echo "🧹 Cleaning up..."
rm -rf "$WORKDIR"

echo "🎉 Job completed at $(date)"