# 1. Point to Library
.libPaths(c("/home/abegunde/R/x86_64-pc-linux-gnu-library/4.5", .libPaths()))

library(qqman)
library(tools)

# 2. CONFIGURATION (ABSOLUTE PATHS)
input_dir  <- "/lscratch/GWAS_practice/DATA/results/summary_statistics"
output_dir <- "/lscratch/GWAS_practice/DATA/results/plots"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# 3. Find files
files <- list.files(input_dir, pattern = "_GWAS_SUMMARY.txt", full.names = TRUE)
print(paste("Found", length(files), "traits to plot."))

# 4. Loop through each file
for (file in files) {
    
    # Get Trait Name
    trait_name <- file_path_sans_ext(basename(file))
    trait_name <- sub("_GWAS_SUMMARY", "", trait_name)
    
    print(paste("Processing:", trait_name))

    # Read Data
    data <- read.table(file, header=TRUE, stringsAsFactors=FALSE)
    data <- na.omit(data)
    
    # Ensure numeric
    data$P <- as.numeric(data$P)
    data$BP <- as.numeric(data$BP)
    data$CHR <- as.numeric(data$CHR)
    
    # Calculate Lambda
    chisq <- qchisq(1 - data$P, 1)
    lambda <- median(chisq) / qchisq(0.5, 1)
    print(paste("  -> Lambda GC:", round(lambda, 4)))

    # --- PLOT 1: MANHATTAN (PDF) ---
    # PDF works on headless servers because it doesn't need pixels.
    pdf_file <- file.path(output_dir, paste0(trait_name, "_Manhattan.pdf"))
    
    pdf(pdf_file, width=12, height=6)
    manhattan(data, 
              main = paste0("Manhattan: ", trait_name, " (L=", round(lambda, 3), ")"), 
              ylim = c(0, 10), 
              col = c("blue4", "orange3"), 
              suggestiveline = -log10(1e-05), 
              genomewideline = -log10(5e-08),
              cex = 0.6)
    dev.off()
    
    # --- PLOT 2: QQ PLOT (PDF) ---
    qq_file <- file.path(output_dir, paste0(trait_name, "_QQ.pdf"))
    
    pdf(qq_file, width=6, height=6)
    qq(data$P, 
       main = paste0("QQ Plot: ", trait_name, " (L=", round(lambda, 3), ")"))
    dev.off()
    
    print(paste("  -> Saved PDF:", trait_name))
}
