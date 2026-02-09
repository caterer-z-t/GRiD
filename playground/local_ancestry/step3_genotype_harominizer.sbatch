#!/bin/bash
#SBATCH --job-name=harmonize
#SBATCH --output=slurm/harmonize.out
#SBATCH --error=slurm/harmonize.err
#SBATCH --time=2-00:00:00
#SBATCH --mem=320G
#SBATCH --partition=general
#SBATCH -n 1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ztcaterer@colorado.edu

# Load modules
module purge
module load java

# Set working directory
WORKDIR=$WORK/genotype_harmonizer
mkdir -p "$WORKDIR"
cd "$WORKDIR" || exit 1

# Define input PLINK paths
VCF2PLINK_DIR="$WORK/vcf2plink"
REFERENCE="$VCF2PLINK_DIR/reference_panel"
TARGET="$VCF2PLINK_DIR/SOL_ch6_LPA_NWD"
GenotypeHarmonizer="$SOFTWARE/GenotypeHarmonizer-1.4.28/GenotypeHarmonizer.jar"

# Run Genotype Harmonizer
java -Xmx30g -jar $GenotypeHarmonizer \
  --input "$TARGET" \
  --inputType PLINK_BED \
  --ref "$REFERENCE" \
  --refType PLINK_BED \
  --output harmonized_SOL_ch6_LPA_NWD \
  --outputType PLINK_BED \
  --update-id

