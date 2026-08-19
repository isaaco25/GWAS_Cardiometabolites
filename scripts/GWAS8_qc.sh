#!/bin/bash
# wakanda_qc_0.2.sh

# ==============================================================================
# 2. CONFIGURATION (DIRECT METHOD)
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

# Structure Cutoff: 0.2 (Requested)
REL_CUTOFF=0.2

source ~/gwas_env/bin/activate
module load plink/1.9 || module load plink

# Wipe folders to prevent "Ghost Files"
rm -rf "$QC_DIR"
mkdir -p "$OUT_DIR" "$QC_DIR"

echo "========================================================"
echo " STARTING FINAL QC V7 (DIRECT GENERATION - 0.2)"
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

# A. Remove PCA Outliers First
# 1. Prune
plink --bfile "$OUT_DIR/temp_samples_clean" --autosome --maf 0.05 --indep-pairwise 50 5 0.2 --out "$QC_DIR/indep_pca" > /dev/null
plink --bfile "$OUT_DIR/temp_samples_clean" --extract "$QC_DIR/indep_pca.prune.in" --make-bed --out "$QC_DIR/pruned_for_pca" > /dev/null
# 2. Calculate PCA
plink --bfile "$QC_DIR/pruned_for_pca" --pca 10 --out "$QC_DIR/pca_check" > /dev/null

# 3. Find Outliers (Python)
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

# B. Apply Relatedness (0.2) - DIRECTLY
# We use the same pruned set from PCA to save time
echo "   -> Applying Relatedness Cutoff (0.2)..."

plink --bfile "$OUT_DIR/temp_samples_no_pca" \
      --extract "$QC_DIR/indep_pca.prune.in" \
      --rel-cutoff "$REL_CUTOFF" \
      --make-bed \
      --out "$QC_DIR/rel_survivors_pruned" > /dev/null

# Now we need to get the LIST of survivors from that pruned file
# and extract them from the FULL dataset.
awk '{print $1, $2}' "$QC_DIR/rel_survivors_pruned.fam" > "$QC_DIR/final_survivor_list.txt"

echo "   -> Survivors after Relatedness: $(wc -l < $QC_DIR/final_survivor_list.txt)"

plink --bfile "$OUT_DIR/temp_samples_no_pca" \
      --keep "$QC_DIR/final_survivor_list.txt" \
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
echo "Original N:     $(wc -l < ${INPUT_DATA}.fam)"
echo "Final N:        $(wc -l < ${OUT_DIR}/wakanda_final_clean.fam)"
echo "Final SNPs:     $(wc -l < ${OUT_DIR}/wakanda_final_clean.bim)"
echo "Cleaned File:   $OUT_DIR/wakanda_final_clean"
echo "========================================================"
