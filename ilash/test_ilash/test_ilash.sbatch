#!/bin/bash
#SBATCH --job-name=test_ilash
#SBATCH --output=slurm/test_ilash.out
#SBATCH --error=slurm/test_ilash.err
#SBATCH --time=11-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=5
#SBATCH --mem=100g
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ztcaterer@colorado.edu

set -euo pipefail
echo "Job started on $(hostname) at $(date)"

# Load modules
module purge
module load gcc/12.2.0
module load ilash/1.0.2
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

# Define paths
WORKDIR=$WORK/test_ilash
OUTPUT_DIR=$LPA/playground/test_ilash

mkdir -p "$WORKDIR" "$OUTPUT_DIR"
cd "$WORKDIR"

# Define Test_File paths
TESTMAP=$SOFTWARE/iLASH/test_files/test.map
TESTPED=$SOFTWARE/iLASH/test_files/test.ped 
CONFIG_FILE=config.txt

# Copy and decompress input
cp "$TESTMAP" .
cp "$TESTPED" .

# Create config
echo "⚙️ Creating iLASH config"
cat > "$CONFIG_FILE" <<EOF
map test.map
ped test.ped
output output.match
slice_size 350
step_size 350
perm_count 20
shingle_size 15
shingle_overlap 0
bucket_count 5
max_thread 2
match_threshold 0.99
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