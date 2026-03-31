import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from scipy.stats import skew

# ==============================================================================
# CONFIGURATION
# ==============================================================================
# File Paths
input_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/harmonized_phenotypes.txt"
output_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/final_qc/phenotypes_clean.txt"
plot_dir = "/home/abegunde/DATA_OMICS_PROJECT/DATA/results/pheno_plots"

# Define your columns based on the headers you provided
ID_COLS = ["FID", "IID"]

# Covariates (We keep these but don't transform them)
COVARIATES = ["Site", "Sex", "Age_wak", "ses_wak", "mvpa_wak", "smoking_status_c_wak", "alcohol_status_c_wak"]

# Traits to Analyze (Quantitative)
# We will check these for outliers and skewness
TRAITS = [
    "standing_height_wak", "weight_wak", "bmi_c_wak", 
    "waist_circumference_wak", "hip_circumference_wak", "waist_hip_r_c_wak",
    "bp_sys_average_wak", "bp_dia_average_wak", 
    "glucose_wak", "insulin_wak", 
    "hdl_wak", "friedewald_ldl_c_c_wak", "cholesterol_1_wak", "triglycerides_wak"
]

# QC Thresholds
SD_THRESHOLD = 4.0  # Remove values > 4 SD from mean (Extreme outliers)
SKEW_THRESHOLD = 1.0 # If skewness > 1.0, apply Log10 transformation

# ==============================================================================
# MAIN SCRIPT
# ==============================================================================
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

print("--- Loading Phenotype Data ---")
# Use 'delim_whitespace=True' because your header looks space/tab separated
df = pd.read_csv(input_file, delim_whitespace=True)

# Ensure numeric columns are actually numeric (coerce errors to NaN)
for col in TRAITS:
    df[col] = pd.to_numeric(df[col], errors='coerce')

print(f"Loaded {len(df)} individuals.")

# ------------------------------------------------------------------------------
# STEP 1: OUTLIER REMOVAL (SD Method)
# ------------------------------------------------------------------------------
print("\n--- Step 1: Removing Extreme Outliers (> 4 SD) ---")
df_clean = df.copy()

for trait in TRAITS:
    mean = df_clean[trait].mean()
    std = df_clean[trait].std()
    
    # Define bounds
    lower = mean - (SD_THRESHOLD * std)
    upper = mean + (SD_THRESHOLD * std)
    
    # Count outliers
    outliers = df_clean[(df_clean[trait] < lower) | (df_clean[trait] > upper)]
    n_out = len(outliers)
    
    if n_out > 0:
        print(f"   -> {trait}: Removing {n_out} outliers (Range: {lower:.2f} to {upper:.2f})")
        # Set outliers to NaN (Missing) instead of dropping the person entirely
        df_clean.loc[(df_clean[trait] < lower) | (df_clean[trait] > upper), trait] = np.nan
    else:
        print(f"   -> {trait}: No outliers found.")

# ------------------------------------------------------------------------------
# STEP 2: NORMALITY CHECK & TRANSFORMATION
# ------------------------------------------------------------------------------
print("\n--- Step 2: Checking Distributions & Transforming ---")
final_traits = []

for trait in TRAITS:
    # Drop NaNs for calculation
    data = df_clean[trait].dropna()
    
    if len(data) == 0:
        print(f"   ! Warning: {trait} has no valid data. Skipping.")
        continue

    # Calculate Skewness
    sk = skew(data)
    
    # Plot Original
    plt.figure(figsize=(10, 4))
    plt.subplot(1, 2, 1)
    sns.histplot(data, kde=True, color="blue")
    plt.title(f"Original: {trait}\nSkew: {sk:.2f}")
    
    # Decide Transformation
    final_col_name = trait
    
    # If highly skewed (e.g., Insulin, Triglycerides), apply LOG
    if abs(sk) > SKEW_THRESHOLD:
        print(f"   -> {trait} is Skewed ({sk:.2f}). Applying LOG10 transformation.")
        
        # Apply Log10 (adding small constant if min <= 0 to avoid errors)
        if df_clean[trait].min() <= 0:
            df_clean[f"LOG_{trait}"] = np.log10(df_clean[trait] + 1)
        else:
            df_clean[f"LOG_{trait}"] = np.log10(df_clean[trait])
            
        final_col_name = f"LOG_{trait}"
        new_sk = skew(df_clean[final_col_name].dropna())
        
        # Plot Transformed
        plt.subplot(1, 2, 2)
        sns.histplot(df_clean[final_col_name], kde=True, color="green")
        plt.title(f"Transformed: LOG_{trait}\nSkew: {new_sk:.2f}")
    else:
        print(f"   -> {trait} is Normal (Skew {sk:.2f}). Keeping original.")
        plt.subplot(1, 2, 2)
        plt.text(0.5, 0.5, "No Transform Needed", ha='center')
        plt.axis('off')

    plt.tight_layout()
    plt.savefig(f"{plot_dir}/hist_{trait}.png")
    plt.close()
    
    final_traits.append(final_col_name)

# ------------------------------------------------------------------------------
# STEP 3: SAVE FINAL FILE
# ------------------------------------------------------------------------------
# Select only ID, Covariates, and Final (Transformed) Traits
cols_to_save = ID_COLS + COVARIATES + [col for col in df_clean.columns if "LOG_" in col or col in TRAITS]

# Remove duplicates if original and LOG versions exist (keep LOG only if transformed)
# Actually, let's keep it simple: Save everything, but the user should use the LOG version if available.
df_clean.to_csv(output_file, sep="\t", index=False, na_rep="NA")

print(f"\n========================================================")
print(f" PHENOTYPE QC COMPLETE")
print(f" Output: {output_file}")
print(f" Plots:  {plot_dir}")
print(f"========================================================")
print("NOTE: For GWAS, use columns starting with 'LOG_' if they exist!")
