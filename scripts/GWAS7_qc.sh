#!/bin/bash
# wakanda_qc_v7.sh

# ==============================================================================
# 2. CONFIGURATION
# ==============================================================================
PROJECT_DIR="/home/abegunde/DATA_OMICS_PROJECT/DATA"
INPUT_DATA="$PROJECT_DIR/cleaned_data/wakanda_harmonized"
OUT_DIR="$PROJECT_DIR/cleaned_data/final_qc"
QC_DIR="$PROJECT_DIR/results/qc_logs"

# Statistical Thresholds
SAMP_MISS_T=0.03
HET_MIN=-0.0106
HET_MAX=0.0152
VAR_MISS_T=0.05
MAF_T=0.01
HWE_STRICT=1e-6

source ~/gwas_env/bin/activate
module load plink/1.9 || module load plink
rm -rf "$QC_DIR"
mkdir -p "$OUT_DIR" "$QC_DIR"

echo "========================================================"
echo " STARTING FINAL QC V7 (DIRECT GENERATION METHOD)"
echo "========================================================"

# ------------------------------------------------------------------------------
# STEP 0: REMOVE DUPLICATE VARIANTS
# ------------------------------------------------------------------------------
echo "--- Step 0: Removing Duplicate Variants ---"
python -c "
import pandas as pd
try:
    df = pd.read_csv('${INPUT_DATA}.bim', sep='\s+', header=None)
    dup_ids = df[df.duplicated(subset=[1], keep='first')]
    dup_ids[1].to_csv('${QC_DIR}/dups.dupvar', index=False, header=False)
    print(f'   -> Found {len(dup_ids)} duplicate variants.')
except Exception as e:
    print(f'   -> Error: {e}')
"

if [ -s "$QC_DIR/dups.dupvar" ]; then
    plink --bfile "$INPUT_DATA" --exclude "$QC_DIR/dups.dupvar" --make-bed --out "$OUT_DIR/temp_nodups" > /dev/null
else
    plink --bfile "$INPUT_DATA" --make-bed --out "$OUT_DIR/temp_nodups" > /dev/null
fi

# ------------------------------------------------------------------------------
# STEP 1: SAMPLE QUALITY
# ------------------------------------------------------------------------------
echo "--- Step 1: Filtering Bad Samples ---"

plink --bfile "$OUT_DIR/temp_nodups" --missing --het --check-sex --out "$QC_DIR/sample_metrics" > /dev/null

# Filter Lists
awk -v min="$HET_MIN" -v max="$HET_MAX" '{if(NR>1 && ($6 < min || $6 > max)) print $1, $2}' "$QC_DIR/sample_metrics.het" > "$QC_DIR/fail_het.txt"
grep "PROBLEM" "$QC_DIR/sample_metrics.sexcheck" | awk '{print $1, $2}' > "$QC_DIR/fail_sex.txt"
awk -v t="$SAMP_MISS_T" '{if(NR>1 && $6 > t) print $1, $2}' "$QC_DIR/sample_metrics.imiss" > "$QC_DIR/fail_miss.txt"

cat "$QC_DIR/fail_het.txt" "$QC_DIR/fail_sex.txt" "$QC_DIR/fail_miss.txt" | sort | uniq > "$QC_DIR/fail_samples_basic.txt"
echo "   -> Removing $(wc -l < $QC_DIR/fail_samples_basic.txt) samples (Het/Sex/Missingness)..."

plink --bfile "$OUT_DIR/temp_nodups" --remove "$QC_DIR/fail_samples_basic.txt" --make-bed --out "$OUT_DIR/temp_samples_clean" > /dev/null

# ------------------------------------------------------------------------------
# STEP 2: STRUCTURE (DIRECT METHOD)
# ------------------------------------------------------------------------------
echo "--- Step 2: Population Structure Checks ---"

# A. PCA Outlier Removal FIRST
# We identify outliers on the clean set, remove them, then run relatedness on the survivors.
plink --bfile "$OUT_DIR/temp_samples_clean" --autosome --maf 0.05 --indep-pairwise 50 5 0.2 --out "$QC_DIR/indep_pca" > /dev/null
plink --bfile "$OUT_DIR/temp_samples_clean" --extract "$QC_DIR/indep_pca.prune.in" --make-bed --out "$QC_DIR/pruned_for_pca" > /dev/null
plink --bfile "$QC_DIR/pruned_for_pca" --pca 10 --out "$QC_DIR/pca_check" > /dev/null

python -c "
import pandas as pd
try:
    df = pd.read_csv('${QC_DIR}/pca_check.eigenvec', sep='\s+', header=None)
    pc1_mean, pc1_sd = df[2].mean(), df[2].std()
    pc2_mean, pc2_sd = df[3].mean(), df[3].std()
    outliers = df[
        (df[2] < pc1_mean - 6*pc1_sd) | (df[2] > pc1_mean + 6*pc1_sd) |
        (df[3] < pc2_mean - 6*pc2_sd) | (df[3] > pc2_mean + 6*pc2_sd)
    ]
    outliers[[0, 1]].to_csv('${QC_DIR}/fail_pca_outliers.txt', sep='\t', index=False, header=False)
    print(f'   -> PCA Outliers Found: {len(outliers)}')
except Exception as e:
    print(f'   -> PCA Error: {e}')
"

echo "   -> Removing PCA Outliers..."
plink --bfile "$OUT_DIR/temp_samples_clean" --remove "$QC_DIR/fail_pca_outliers.txt" --make-bed --out "$OUT_DIR/temp_samples_no_pca" > /dev/null

# B. Relatedness (0.9) - DIRECT GENERATION
# We run rel-cutoff AND make-bed in the same command. This bypasses list file errors.
# We also create a new pruned set for this specific file to be safe.
echo "   -> Running Relatedness Pruning (0.9)..."

plink --bfile "$OUT_DIR/temp_samples_no_pca" --extract "$QC_DIR/indep_pca.prune.in" --make-bed --out "$QC_DIR/pruned_for_rel" > /dev/null

# CRITICAL: This mimics the plot script exactly. 
# It reads 'pruned_for_rel', calculates relatedness, applies cutoff 0.9, 
# and writes the CLEAN list of survivors to 'wakanda_samples_final'.
plink --bfile "$QC_DIR/pruned_for_rel" \
      --rel-cutoff 0.9 \
      --make-bed \
      --out "$OUT_DIR/wakanda_samples_final_pruned" > /dev/null

# Note: The output above is pruned (few SNPs). We need the FULL SNPs for the final file.
# So we extract the list of SURVIVORS from the .fam file we just made.
awk '{print $1, $2}' "$OUT_DIR/wakanda_samples_final_pruned.fam" > "$QC_DIR/final_survivors.txt"

# Now extract those survivors from the full dataset
plink --bfile "$OUT_DIR/temp_samples_no_pca" \
      --keep "$QC_DIR/final_survivors.txt" \
      --make-bed \
      --out "$OUT_DIR/wakanda_samples_final" > /dev/null

# ------------------------------------------------------------------------------
# STEP 3: VARIANT QUALITY
# ------------------------------------------------------------------------------
echo "--- Step 3: Final Variant Filters ---"

plink --bfile "$OUT_DIR/wakanda_samples_final" \
      --autosome \
      --geno "$VAR_MISS_T" \
      --maf "$MAF_T" \
      --hwe "$HWE_STRICT" \
      --make-bed \
      --out "$OUT_DIR/wakanda_final_clean"

echo "========================================================"
echo " QC SUCCESSFUL"
echo "========================================================"
echo "Final Data:     $(wc -l < ${OUT_DIR}/wakanda_final_clean.fam) individuals"
echo "Cleaned File:   $OUT_DIR/wakanda_final_clean"
echo "========================================================"
