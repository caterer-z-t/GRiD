#!/bin/bash
module load anaconda/4.3.0
source /nas/longleaf/apps/anaconda/4.3.0/anaconda/etc/profile.d/conda.sh
conda activate lpa

# You can add any commands to run here
echo "✅ Conda environment activated!"
python --version
