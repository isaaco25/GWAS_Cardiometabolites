import pandas as pd

# --- CONFIGURATION ---
#  files names to be replaced with their correct path 
fam_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/wakanda.fam"   #  PLINK .fam file
pheno_file = "/home/abegunde/DATA_OMICS_PROJECT/DATA/Phenotypes.txt"    # phenotype .txt file
output_keep_list = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/samples_to_keep.txt" # samples to be kept
output_clean_pheno = "/home/abegunde/DATA_OMICS_PROJECT/DATA/cleaned_data/harmonized_phenotypes.txt" # new phenotypic harmonised file

# --- 1. LOAD GENOTYPE IDs ---
# .fam files are space-delimited and have no header.
# We only need the first two columns: FID (Family ID) and IID (Individual ID).
geno_df = pd.read_csv(fam_file, delim_whitespace=True, header=None, usecols=[0, 1], names=["FID", "IID"], dtype=str)

# --- 2. LOAD PHENOTYPE DATA ---
# Adjust 'sep' if your file is comma-separated (sep=',')
# We assume your text file has headers like "FID", "IID", "LDL", etc.
pheno_df = pd.read_csv(pheno_file, sep=r'\s+', dtype=str)

# --- 3. HARMONIZE (FIND INTERSECTION) 
# This merges the two dataframes based on FID and IID, keeping only matching rows.
common_df = pd.merge(geno_df, pheno_df, on=["FID", "IID"], how="inner")

# --- 4. EXPORT FILES ---

# A. Save the Keep List for PLINK (Only FID and IID required)
common_df[["FID", "IID"]].to_csv(output_keep_list, sep="\t", index=False, header=False)

# B. Save the Harmonized Phenotype File (Clean data for analysis)
common_df.to_csv(output_clean_pheno, sep="\t", index=False)

# --- 5. REPORT ---
print(f"Total Genotyped Samples: {len(geno_df)}")
print(f"Total Phenotyped Samples: {len(pheno_df)}")
print(f"------------------------------------------------")
print(f"Intersection (Samples in both): {len(common_df)}")

if len(common_df) == 0:
    print("WARNING: No matches found!.")
else:
    print(f"Files '{output_keep_list}' and '{output_clean_pheno}' created successfully.")
