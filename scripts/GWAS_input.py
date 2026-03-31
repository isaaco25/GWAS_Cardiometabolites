import pandas as pd
import os

# ==============================================================================
# CONFIGURATION
# ==============================================================================
pheno_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/final_qc/phenotypes_clean.txt"
pca_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/results/population_structure/final_pca.eigenvec"
out_dir = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/GWAS_ready"

if not os.path.exists(out_dir):
    os.makedirs(out_dir)

print("--- Preparing GWAS Input Files ---")

# ------------------------------------------------------------------------------
# 1. LOAD & CLEAN PHENOTYPES
# ------------------------------------------------------------------------------
df_pheno = pd.read_csv(pheno_file, sep="\t")
print(f"Loaded {len(df_pheno)} phenotype records.")

# IDENTIFY & DROP OLD PCs
# We find any column that starts with "PC" (e.g., PC1, PC2, PC1_x, PC_old)
old_pcs = [col for col in df_pheno.columns if col.startswith("PC")]

if old_pcs:
    print(f"-> Found {len(old_pcs)} existing PC columns (e.g., {old_pcs[:3]}...). Deleting them.")
    df_pheno = df_pheno.drop(columns=old_pcs)
else:
    print("-> No existing PC columns found in phenotype file.")

# ------------------------------------------------------------------------------
# 2. LOAD NEW PCA
# ------------------------------------------------------------------------------
# PLINK .eigenvec has NO header. We must assign names manually.
try:
    df_pca = pd.read_csv(pca_file, sep="\s+", header=None)
    
    # Check shape to ensure we have FID, IID, + 10 PCs
    num_cols = df_pca.shape[1]
    if num_cols >= 12:
        # Standard: FID, IID, PC1...PC10...
        # We only keep top 10 even if file has 20
        df_pca = df_pca.iloc[:, 0:12]
        df_pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, 11)]
    else:
        # Fallback for weird shapes
        print(f"WARNING: PCA file has {num_cols} columns. Using available columns.")
        df_pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, num_cols-1)]
    
    print(f"Loaded {len(df_pca)} PCA records from PLINK.")
    
except Exception as e:
    print(f"CRITICAL ERROR loading PCA file: {e}")
    exit(1)

# ------------------------------------------------------------------------------
# 3. MERGE
# ------------------------------------------------------------------------------
# Inner join: Keep only individuals present in BOTH files
df_merged = pd.merge(df_pheno, df_pca, on=["FID", "IID"])
print(f"Merged Dataset: {len(df_merged)} individuals with Genotype + Phenotype.")

# ------------------------------------------------------------------------------
# 4. DEFINE COLUMNS
# ------------------------------------------------------------------------------
# Covariates: Fixed + The NEW Genetic PCs
covar_cols = ["FID", "IID", "Sex", "Age_wak", "Site", "ses_wak", "mvpa_wak"] + [f"PC{i}" for i in range(1, 11)]

# Phenotypes: LOG columns + Original Quantitative traits
# We exclude the raw versions of the transformed traits to avoid confusion
all_cols = df_merged.columns.tolist()

pheno_cols = ["FID", "IID"] + [
    c for c in all_cols 
    if (c not in covar_cols) and ("LOG_" in c or c in [
        "bmi_c_wak", "standing_height_wak", "weight_wak", 
        "waist_circumference_wak", "hip_circumference_wak", 
        "bp_sys_average_wak", "bp_dia_average_wak", "glucose_wak", 
        "friedewald_ldl_c_c_wak", "cholesterol_1_wak"
    ])
]

# ------------------------------------------------------------------------------
# 5. SAVE FILES
# ------------------------------------------------------------------------------
try:
    # Save Covariates
    df_merged[covar_cols].to_csv(f"{out_dir}/gwas_covariates.txt", sep="\t", index=False, na_rep="NA")
    print(f"-> Created Covariate File: {out_dir}/gwas_covariates.txt")

    # Save Phenotypes
    df_merged[pheno_cols].to_csv(f"{out_dir}/gwas_phenotypes.txt", sep="\t", index=False, na_rep="NA")
    print(f"-> Created Phenotype File: {out_dir}/gwas_phenotypes.txt")
    
except KeyError as e:
    print(f"\nERROR: Columns missing during save. \nExpected: {covar_cols}\nFound in data: {df_merged.columns.tolist()}")
    exit(1)

print("\n[READY] Traits available for analysis:")
for c in pheno_cols[2:]:
    print(f" - {c}")
