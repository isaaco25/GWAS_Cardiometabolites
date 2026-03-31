#!/bin/bash
#SBATCH --job-name=Final_GWAS
#SBATCH --mem=4gb
#SBATCH --cpus-per-task=1
#SBATCH --array=1-22           # Run one job per chromosome
#SBATCH --output=logs/FinalGWAS_%a.log

# ==========================================
# CONFIGURATION
# ==========================================

# 1. INPUT: The merged PLINK binary file 

BFILE="cleaned_data/imputation_ready_strict_fix/final_merged_data/wakanda_imputed_final"

# 2. PHENOTYPES: Path to phenotype text file
PHENO_FILE="cleaned_data/GWAS_ready/gwas_phenotypes.txt"

# 3. COVARIATES: Path to covariate file
COVAR_FILE="cleaned_data/GWAS_ready/gwas_covariates.txt"

# 4. THE MODEL: Exact column names from text file
COVAR_NAMES="Age_wak,Sex,PC1,PC2,PC3"

# 5. OUTPUT: save the results
OUT_DIR="results/final_gwas_results"

module load plink/1.9 || module load plink

mkdir -p "$OUT_DIR"

# =====================================
# SAFETY CHECKS
# =====================================
if [ ! -f "${BFILE}.bed" ]; then
    echo "CRITICAL ERROR: Input file not found at ${BFILE}.bed"
    exit 1
fi
if [ ! -f "$PHENO_FILE" ]; then
    echo "CRITICAL ERROR: Phenotype file not found at $PHENO_FILE"
    exit 1
fi
if [ ! -f "$COVAR_FILE" ]; then
    echo "CRITICAL ERROR: Covariate file not found at $COVAR_FILE"
    exit 1
fi

# ====================================
# THE GWAS COMMAND
# ====================================
CHR=$SLURM_ARRAY_TASK_ID
echo "Running Final GWAS for Chromosome $CHR..."

# EXPLANATION OF FLAGS:
# --bfile:        Input data
# --chr:          Only analyze this specific chromosome (1-22)
# --pheno:        Load the phenotype file
# --all-pheno:    Run a separate GWAS for EVERY column in the pheno file
# --covar:        Load the covariate file
# --covar-name:   Only use Age, Sex, PC1-3 (ignore PC4-10)
# --linear hide-covar: Use Linear Regression but don't output lines for Age/Sex
# --adjust:       Calculate FDR and Bonferroni corrections
# --ci 0.95:      Calculate 95% Confidence Intervals

plink \
    --threads $SLURM_CPUS_PER_TASK \
    --bfile "$BFILE" \
    --chr "$CHR" \
    --pheno "$PHENO_FILE" \
    --all-pheno \
    --covar "$COVAR_FILE" \
    --covar-name "$COVAR_NAMES" \
    --linear hide-covar \
    --ci 0.95 \
    --allow-no-sex \
    --adjust \
    --memory 15000 \
    --out "${OUT_DIR}/chr${CHR}_results"

echo "Chromosome $CHR Complete."
