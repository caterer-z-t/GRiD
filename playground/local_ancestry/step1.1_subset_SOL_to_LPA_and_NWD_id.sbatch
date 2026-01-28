#!/bin/bash
#SBATCH --job-name=subset_sol
#SBATCH --output=slurm/%x.out
#SBATCH --error=slurm/%x.err
#SBATCH --time=2-00:00:00
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mem=100G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load required modules
module purge
module load plink/2.00

# Subset SOL data to LPA \pm 50mb

# Define file paths
SOL_DIR="${WORK}/vcf2plink"

if [ ! -f "${SOL_DIR}/SOL_ch6_LPA.bim" ]; then
    plink2 --bfile "${SOL_DIR}/SOL_ch6" \
        --chr 6 \
        --from-bp 135000000 \
        --to-bp 185000000 \
        --make-bed \
        --out "${SOL_DIR}/SOL_ch6_LPA"
fi

# Subset to NWD and LPA IDs
# Get a list of NWD_IDs from $CRAM directory
for f in ${CRAM}/NWD*.cram; do 
    basename "${f}" .b38.irc.v1_subset.cram >> "${SOL_DIR}/NWD_IDs.txt"
done

# Make PLINK keep file (two columns)
awk '{print $1, $1}' "${SOL_DIR}/NWD_IDs.txt" >> "${SOL_DIR}/NWD_IDs_keep.txt"

if [ ! -f "${SOL_DIR}/SOL_ch6_LPA_NWD.bim" ]; then
    plink2 --bfile "${SOL_DIR}/SOL_ch6_LPA" \
        --keep "${SOL_DIR}/NWD_IDs_keep.txt" \
        --make-bed \
        --out "${SOL_DIR}/SOL_ch6_LPA_NWD"
fi