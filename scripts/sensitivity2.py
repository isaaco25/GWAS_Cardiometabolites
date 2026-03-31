import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# CONFIG
dir_path = "/home/abegunde/DATA_OMICS_PROJECT/DATA/results/sensitivity"
sns.set_style("whitegrid")

# Helper function for consistent plotting
def create_sensitivity_plot(x, y, xlabel, ylabel, title, filename, vlines=None, log_x=False):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y, marker='o', linestyle='-', color='teal', linewidth=2)
    
    if log_x:
        plt.xscale('log')
        
    plt.title(title, fontsize=14)
    plt.xlabel(xlabel, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    
    # Add vertical guide lines if requested (e.g., standard thresholds)
    if vlines:
        for val, label, color in vlines:
            plt.axvline(x=val, color=color, linestyle='--', label=label)
        plt.legend()
        
    plt.grid(True, which="both", linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot: {filename}")

# ==========================================
# 1. RELATEDNESS PLOT
# ==========================================
try:
    rel_file = f"{dir_path}/sensitivity_relatedness.txt"
    if os.path.exists(rel_file):
        df_rel = pd.read_csv(rel_file, sep="\t")
        
        # Guide lines for standard degrees of relationship
        guides = [
            (0.125, "3rd Degree (Cousins)", "red"),
            (0.2, "2nd Degree (Half-Sibs)", "orange")
        ]
        
        create_sensitivity_plot(
            df_rel['Threshold'], df_rel['Survivors'],
            "Relatedness Threshold (PI_HAT)", "Survivors (N)",
            "Sample Retention vs Relatedness Strictness",
            "plot_relatedness_curve.png", guides
        )
    else:
        print(f"Skipping Relatedness: {rel_file} not found.")
except Exception as e:
    print(f"Error plotting Relatedness: {e}")

# ==========================================
# 2. MAF PLOT (Minor Allele Frequency)
# ==========================================
try:
    maf_file = f"{dir_path}/raw_qc_metrics.frq"
    if os.path.exists(maf_file):
        # Read .frq file (whitespace delimited)
        df_maf = pd.read_csv(maf_file, delim_whitespace=True)
        
        # Test specific thresholds
        thresholds = [0.001, 0.005, 0.01, 0.02, 0.03, 0.05]
        survivors = [len(df_maf[df_maf['MAF'] >= t]) for t in thresholds]
        
        guides = [
            (0.01, "1% Cutoff (Rare)", "red"),
            (0.05, "5% Cutoff (Common)", "blue")
        ]
        
        create_sensitivity_plot(
            thresholds, survivors,
            "MAF Threshold", "SNPs Remaining",
            "Variant Retention vs MAF Threshold",
            "plot_maf_curve.png", guides
        )
    else:
        print(f"Skipping MAF: {maf_file} not found.")
except Exception as e:
    print(f"Error plotting MAF: {e}")

# ==========================================
# 3. HWE PLOT (Hardy-Weinberg Equilibrium)
# ==========================================
try:
    hwe_file = f"{dir_path}/raw_qc_metrics.hwe"
    if os.path.exists(hwe_file):
        df_hwe = pd.read_csv(hwe_file, delim_whitespace=True)
        
        # HWE logic: We keep SNPs with P > Threshold
        # We test log-scale thresholds
        thresholds = [1e-20, 1e-10, 1e-6, 1e-5, 1e-4, 1e-3]
        survivors = [len(df_hwe[df_hwe['P'] > t]) for t in thresholds]
        
        guides = [(1e-6, "Standard (1e-6)", "red")]
        
        create_sensitivity_plot(
            thresholds, survivors,
            "HWE P-value Threshold (Log Scale)", "SNPs Remaining",
            "Variant Retention vs HWE Strictness",
            "plot_hwe_curve.png", guides, log_x=True
        )
    else:
        print(f"Skipping HWE: {hwe_file} not found.")
except Exception as e:
    print(f"Error plotting HWE: {e}")
