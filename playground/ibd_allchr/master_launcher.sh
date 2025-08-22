#!/bin/bash

CHROM=$1

if [[ -z "$CHROM" ]]; then
  echo "Usage: $0 <chromosome_number>"
  exit 1
fi

export CHROM=$CHROM

# Step 1: VCF to PLINK
jid1=$(sbatch --export=CHROM=$CHROM step1_vcf_to_plink.sh | awk '{print $4}')
echo "Step 1 submitted: Job ID = $jid1"

# Step 2: Eagle phasing
jid2=$(sbatch --dependency=afterok:$jid1 --export=CHROM=$CHROM step2_eagle.sh | awk '{print $4}')
echo "Step 2 submitted: Job ID = $jid2"

# Step 3: iLASH
jid3=$(sbatch --dependency=afterok:$jid2 --export=CHROM=$CHROM step3_ilash.sh | awk '{print $4}')
echo "Step 3 submitted: Job ID = $jid3"
