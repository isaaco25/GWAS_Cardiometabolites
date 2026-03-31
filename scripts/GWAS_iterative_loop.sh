#!/bin/bash
# run_iterative_gwas.sh
#SBATCH --job-name=isaac_abegunde   # Job name
#SBATCH --mail-type=BEGIN,END,FAIL  # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=abegunde.isaaco@pg.funaab.edu.ng
#SBATCH --mem=64gb                     # Job memory request
#SBATCH --cpus-per-task=4             # Request 4 CPU cores
#SBATCH --ntasks=4
#SBATCH --nodelist=node01
#SBATCH --output=stepwise_gwas.log      # Standard output log
#SBATCH --error=stepwise_gwas.err       # Error log
# ==============================================================================
# CONFIGURATION
# ==============================================================================
PROJECT_DIR="/lscratch/GWAS_practice/DATA"
GENO_DATA="$PROJECT_DIR/cleaned_data/final_qc/wakanda_final_clean"
PHENO_FILE="$PROJECT_DIR/cleaned_data/GWAS_ready/gwas_phenotypes.txt"
COVAR_FILE="$PROJECT_DIR/cleaned_data/GWAS_ready/gwas_covariates.txt"
OUT_DIR="$PROJECT_DIR/results/gwas_results/stepwise_analysis"

mkdir -p "$OUT_DIR"

module load plink/1.90
# DEFINE YOUR ORDERED COVARIATE LIST
# The script will add these one by one in this exact order.
# We usually start with demographics, then SES/Lifestyle, then Genetic PCs.
COVAR_ORDER=("Sex" "Age_wak" "PC1" "PC2" "PC3" "PC4" "PC5" "PC6" "PC7" "PC8" "PC9" "PC10")

# DETECT TRAITS (Same as before)
TRAITS=$(python -c "import pandas as pd; df=pd.read_csv('$PHENO_FILE', sep='\t', nrows=0); print(' '.join(df.columns[2:]))")

echo "========================================================"
echo " STARTING STEPWISE CUMULATIVE GWAS"
echo "========================================================"

# Initialize empty string for the cumulative list
CURRENT_COVAR_STRING=""
STEP_COUNT=0

# ------------------------------------------------------------------------------
# OUTER LOOP: Iterate through covariates one by one
# ------------------------------------------------------------------------------
for NEW_COVAR in "${COVAR_ORDER[@]}"; do
    STEP_COUNT=$((STEP_COUNT+1))
    
    # Add the new covariate to the comma-separated string
    if [ -z "$CURRENT_COVAR_STRING" ]; then
        CURRENT_COVAR_STRING="$NEW_COVAR"
    else
        CURRENT_COVAR_STRING="$CURRENT_COVAR_STRING,$NEW_COVAR"
    fi

    echo "--------------------------------------------------------"
    echo " STEP $STEP_COUNT: Model = Phenotype ~ Genotype + [ $CURRENT_COVAR_STRING ]"
    echo "--------------------------------------------------------"

    # Create a sub-folder for this step to keep things organized
    STEP_DIR="$OUT_DIR/step${STEP_COUNT}_added_${NEW_COVAR}"
    mkdir -p "$STEP_DIR"

    # --------------------------------------------------------------------------
    # INNER LOOP: Run GWAS for every trait using the current covariate list
    # --------------------------------------------------------------------------
    for TRAIT in $TRAITS; do
        # echo "   -> Analyzing $TRAIT..."
        
        plink --bfile "$GENO_DATA" \
              --pheno "$PHENO_FILE" \
              --threads 4 \
              --pheno-name "$TRAIT" \
              --covar "$COVAR_FILE" \
              --covar-name "$CURRENT_COVAR_STRING" \
              --linear hide-covar \
              --adjust \
              --allow-no-sex \
              --out "$STEP_DIR/assoc_$TRAIT" > /dev/null
              
        # I added '> /dev/null' to keep the terminal clean, remove it if you want to see PLINK logs
    done
    
    echo "   -> Completed all traits for Step $STEP_COUNT."
done

echo "========================================================"
echo " STEPWISE ANALYSIS COMPLETE."
echo " Check folder: $OUT_DIR"
echo "========================================================"
