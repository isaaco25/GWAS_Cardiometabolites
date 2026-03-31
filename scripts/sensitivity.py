import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# CONFIG
out_dir = "/home/abegunde/DATA_OMICS_PROJECT/DATA/results/sensitivity"
prefix = f"{out_dir}/raw_qc_metrics"

# Setup plots style
sns.set_style("whitegrid")

def calculate_survivors(df, col, thresholds, less_than=True):
    results = []
    total = len(df)
    for t in thresholds:
        if less_than:
            count = len(df[df[col] < t])
        else:
            count = len(df[df[col] > t])
        results.append({'Threshold': t, 'Survivors': count, 'Percent': (count/total)*100})
    return pd.DataFrame(results)

# ==========================================
# 1. SAMPLE MISSINGNESS (F_MISS)
# ==========================================
try:
    print("Analyzing Sample Missingness...")
    imiss = pd.read_csv(f"{prefix}.imiss", delim_whitespace=True)
    
    # Define thresholds to check
    thresh_imiss = [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.07, 0.10]
    df_imiss_sens = calculate_survivors(imiss, "F_MISS", thresh_imiss, less_than=True)
    
    # Save Table
    df_imiss_sens.to_csv(f"{out_dir}/table_sample_missingness.txt", sep="\t", index=False)
    
    # Plot Distribution + Threshold Curve
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    # Histogram
    sns.histplot(imiss['F_MISS'], bins=50, ax=ax[0], color='skyblue')
    ax[0].set_title("Distribution of Sample Missingness")
    ax[0].set_xlabel("Missingness Fraction")
    
    # Threshold Curve
    ax[1].plot(df_imiss_sens['Threshold'], df_imiss_sens['Survivors'], marker='o', color='navy')
    ax[1].set_title("Survivors vs Threshold")
    ax[1].set_ylabel("Number of Samples Kept")
    ax[1].set_xlabel("Missingness Threshold (<)")
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/plot_sample_missingness.png")
    plt.close()
except Exception as e:
    print(f"Skipped Sample Miss: {e}")

# ==========================================
# 2. VARIANT MISSINGNESS (F_MISS)
# ==========================================
try:
    print("Analyzing Variant Missingness...")
    lmiss = pd.read_csv(f"{prefix}.lmiss", delim_whitespace=True)
    
    thresh_lmiss = [0.01, 0.02, 0.03, 0.05, 0.10]
    df_lmiss_sens = calculate_survivors(lmiss, "F_MISS", thresh_lmiss, less_than=True)
    
    df_lmiss_sens.to_csv(f"{out_dir}/table_variant_missingness.txt", sep="\t", index=False)
    
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    sns.histplot(lmiss['F_MISS'], bins=50, ax=ax[0], color='salmon', log_scale=(False, True)) # Log scale y for variants often helps
    ax[0].set_title("Distribution of SNP Missingness (Log Y)")
    
    ax[1].plot(df_lmiss_sens['Threshold'], df_lmiss_sens['Survivors'], marker='o', color='darkred')
    ax[1].set_title("SNPs vs Threshold")
    
    plt.tight_layout()
    plt.savefig(f"{out_dir}/plot_variant_missingness.png")
    plt.close()
except Exception as e:
    print(f"Skipped Variant Miss: {e}")

# ==========================================
# 3. HETEROZYGOSITY (SD Approach)
# ==========================================
try:
    print("Analyzing Heterozygosity...")
    het = pd.read_csv(f"{prefix}.het", delim_whitespace=True)
    
    mean_f = het['F'].mean()
    sd_f = het['F'].std()
    
    # Check SD thresholds (keeping samples WITHIN +/- X SDs)
    sd_thresholds = [2, 2.5, 3, 3.5, 4, 4.5, 5]
    het_results = []
    
    for sd in sd_thresholds:
        lower = mean_f - (sd * sd_f)
        upper = mean_f + (sd * sd_f)
        count = len(het[(het['F'] >= lower) & (het['F'] <= upper)])
        het_results.append({'SD_Cutoff': sd, 'Survivors': count, 'Lower_F': lower, 'Upper_F': upper})
        
    df_het_sens = pd.DataFrame(het_results)
    df_het_sens.to_csv(f"{out_dir}/table_heterozygosity.txt", sep="\t", index=False)
    
    # Plot
    plt.figure(figsize=(8, 6))
    sns.histplot(het['F'], bins=50, color='purple', kde=True)
    plt.axvline(mean_f - (3*sd_f), color='r', linestyle='--', label='Mean +/- 3SD')
    plt.axvline(mean_f + (3*sd_f), color='r', linestyle='--')
    plt.title(f"Heterozygosity Distribution (Mean: {mean_f:.4f}, SD: {sd_f:.4f})")
    plt.legend()
    plt.savefig(f"{out_dir}/plot_heterozygosity.png")
    plt.close()
except Exception as e:
    print(f"Skipped Het: {e}")

print("Done! Check the 'results/sensitivity' folder for plots and tables.")
