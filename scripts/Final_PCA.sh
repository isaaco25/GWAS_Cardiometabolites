#!/bin/bash
# wakanda_final_pca.sh

# ==============================================================================
# CONFIGURATION
# ==============================================================================
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/final_qc/wakanda_final_clean"
OUT_DIR="$PROJECT_DIR/results/population_structure"

mkdir -p "$OUT_DIR"
module load plink/1.9 || module load plink


echo "========================================================"
echo " GENERATING FINAL PCA (COVARIATES)"
echo "========================================================"

# 1. Prune the clean dataset (LD Pruning)
# We prune again because PCA requires independent markers
plink --bfile "$INPUT_DATA" \
      --indep-pairwise 50 5 0.2 \
      --out "$OUT_DIR/clean_indep"

# 2. Calculate PCA (Top 10 PCs)
plink --bfile "$INPUT_DATA" \
      --extract "$OUT_DIR/clean_indep.prune.in" \
      --pca 10 \
      --out "$OUT_DIR/final_pca"

echo "PCA Calculated. Eigenvectors saved to: $OUT_DIR/final_pca.eigenvec"

# 3. Visualization (Python)
# We plot PC1 vs PC2 to visualize your final cohort structure
source ~/gwas_env/bin/activate

python -c "
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

pca_file = '${OUT_DIR}/final_pca.eigenvec'
out_plot = '${OUT_DIR}/plot_final_pca.png'

try:
    # Load Data (PLINK .eigenvec has no header)
    # Col 0=FID, 1=IID, 2=PC1, 3=PC2 ...
    df = pd.read_csv(pca_file, sep='\s+', header=None)
    
    plt.figure(figsize=(10, 8))
    sns.scatterplot(x=df[2], y=df[3], alpha=0.6, color='darkblue', edgecolor=None)
    
    plt.title('Final Population Structure (PC1 vs PC2)', fontsize=14)
    plt.xlabel(f'Principal Component 1', fontsize=12)
    plt.ylabel(f'Principal Component 2', fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()
    
    plt.savefig(out_plot)
    print(f'PCA Plot saved: {out_plot}')
    
except Exception as e:
    print(f'Error plotting PCA: {e}')
"

echo "========================================================"
echo " PROCESS COMPLETE"
echo " You are now ready for Association Testing."
echo "========================================================"
