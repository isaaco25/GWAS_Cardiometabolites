#!/bin/bash
#SBATCH --job-name=isaac_abegunde_process
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2            # PLINK 2 is fast, 2-4 cpus is usually plenty
#SBATCH --mem=8G                     # 8GB is usually safe for per-chr jobs
#SBATCH --array=1-22                 # <--- This runs the job 22 times (Chr 1 to 22)
#SBATCH --output=logs/process_michigan.log  # Saves nice log files per chromosome


# The chromosome number for this specific task
CHR=${SLURM_ARRAY_TASK_ID}

# Where are your input VCFs located?
IN_DIR="cleaned_data/imputation_ready_strict_fix/michigan_result"

# Where should the clean output go?
OUT_DIR="cleaned_data/imputation_ready_strict_fix/michigan_result/imputation_clean"

# Create output directory if it doesn't exist
mkdir -p ${OUT_DIR}
mkdir -p logs


module purge
# Intentionally plink2, not plink/1.9 — needed for native VCF/dosage import post-imputation
module load plink/2.0  

# Verify we have the right version in the log file
echo "Running on node: $(hostname)"
echo "PLINK version check:"
plink2 --version

# ------------
# 3. EXECUTION
# -------------

INPUT_FILE="${IN_DIR}/chr${CHR}.dose.vcf.gz"
OUTPUT_PREFIX="${OUT_DIR}/chr${CHR}_imputed_clean"

# SAFETY CHECK: Does the input file actually exist?
if [ ! -f "$INPUT_FILE" ]; then
    echo "ERROR: Input file not found: $INPUT_FILE"
    exit 1
fi

echo "Processing Chromosome ${CHR}..."

# Run PLINK 2
# Note: I added '1>' and '2>' to capture explicit logs into the log folder
plink2 \
    --vcf ${INPUT_FILE} \
    --double-id \
    --allow-extra-chr \
    --extract-if-info R2 '>' 0.3 \
    --maf 0.01 \
    --new-id-max-allele-len 50 \
    --vcf-half-call m \
    --make-bed \
    --out ${OUTPUT_PREFIX} \
    --memory 7000  # Request slightly less than the SBATCH --mem limit

# ----------------
# 4. VERIFICATION
# ---------------

if [ -f "${OUTPUT_PREFIX}.bed" ]; then
    echo "SUCCESS: Chromosome ${CHR} finished correctly."
else
    echo "FAILURE: The .bed file was not created for Chromosome ${CHR}."
    exit 1
fi
