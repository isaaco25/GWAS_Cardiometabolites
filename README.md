# GWAS_Cardiometabolites
# GWAS Analysis of Cardiometabolic Phenotypes

Genome-Wide Association Study (GWAS) pipeline for identifying genetic variants associated with cardiometabolic traits in a human population cohort.

## Project Overview

This project was developed as part of the **International Certificate in Bioinformatics and Genomics (CIBiG)** training program at the WAVE Centre, Abidjan (November–December 2025).

### Research Question
Identify genetic variants associated with cardiometabolic traits (BMI, blood pressure, hypertension) in a population cohort using a GWAS framework.

### Study Design
- **Sample Size:** 2,063 individuals (1,416 post-QC harmonization)
- **Genotype Format:** PLINK (.bed, .bim, .fam)
- **Reference Genome:** GRCh37
- **Genotyping Platform:** H3Africa array (Illumina)
- **Phenotypes:** Anthropometric, metabolic, cardiovascular, lifestyle traits

## Repository Structure
```
├── scripts/           # Analysis scripts (Bash, R)
└── README.md
```

## Analytical Workflow

0. **Setup** — Create the `logs/` directory (`mkdir -p logs`) before submitting any SLURM script 
1. **Sample Harmonization** — Match genotype and phenotype data
2. **Genotype QC (PLINK)** — Filter by missingness, MAF, HWE
3. **Population Structure** — PCA for ancestry correction (PC1–PC10)
4. **VCF Generation** — Convert autosomes to bgzipped VCFs
5. **Chromosome X Processing** — PAR splitting, sex checks
6. **Genotype Imputation** — Michigan Imputation Server
7. **Post-Imputation QC** — INFO/R² threshold, MAF filtering
8. **GWAS Analysis** — Linear regression (additive model)
9. **Visualization** — Manhattan plots, Q-Q plots

## Tools & Software

| Tool | Version | Purpose |
|------|---------|---------|
| PLINK | 1.9 / 2.0 | QC, PCA, association testing |
| bcftools | - | VCF file manipulation |
| tabix | - | VCF indexing |
| vcftools | - | VCF filtering and statistics |
| GATK | - | Variant processing |
| FastQC | - | Sequence quality control |
| Trimmomatic | - | Read trimming |
| Michigan Imputation Server | - | Genotype imputation |
| FUMA-GWAS | - | Functional annotation |
| R | 4.5.2 | Statistical analysis & visualization |

## Computational Environment

- **HPC Cluster:** WAVE Centre / SLURM workload manager
- **Environment Management:** Conda

## Author

**Isaac Abegunde** — [GitHub](https://github.com/isaaco25)

## Acknowledgments

- CIBiG Training Program
- WAVE Centre, Abidjan
