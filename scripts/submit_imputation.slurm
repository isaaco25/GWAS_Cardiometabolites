#!/bin/bash
#SBATCH --job-name=Imp_Upload
#SBATCH --mail-type=ALL
#SBATCH --output=imputation_upload.log
#SBATCH --error=imputation_upload.err
#SBATCH --mem=8gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2

# 1. Run the Imputation Bot
/home/abegunde/imputationbot impute \
  --files cleaned_data/imputation_ready_strict_fix/*.vcf.gz \
  --refpanel 1000g-phase-3-v5 \
  --population afr \
  --name "Isaac_Imputation_Run" \


