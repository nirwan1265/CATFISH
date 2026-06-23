# ==============================================================================
# Arabidopsis Total Nitrogen (TN) 0-5cm Mean - SNP to Gene Annotation
# Phenotype: nitrogen_0_5cm_mean_1000
# ==============================================================================
#
# This script:
#   1. Loads processed GWAS results (arabidopsis_TN_0_5_mean.csv)
#   2. Maps SNPs to genes using 25kb upstream/downstream window
#   3. Outputs annotated gene-level statistics
#
# ==============================================================================

library(data.table)
library(dplyr)
library(ggplot2)

if (!requireNamespace("qqman", quietly = TRUE)) {
  install.packages("qqman")
}
library(qqman)

# ==============================================================================
# USER PARAMETERS
# ==============================================================================

WINDOW_SIZE <- 25000  # 25kb upstream and downstream

# Input files
GWAS_FILE <- "/Users/nirwantandukar/Documents/Research/data/Arabidopsis/raw_gwas/arabidopsis_TN_0_5_mean.csv"
GENE_LOC_FILE <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/at.genes.loc"

# Output directory
OUT_DIR <- "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/GWAS"

# Create output directory if it doesn't exist
if (!dir.exists(OUT_DIR)) {
  dir.create(OUT_DIR, recursive = TRUE)
}

# ==============================================================================
# 1. Load Gene Location File
# ==============================================================================

cat("Loading gene location file...\n")
gene_loc <- fread(GENE_LOC_FILE)
cat("Total genes in annotation:", nrow(gene_loc), "\n")

# Create extended gene regions (25kb upstream and downstream)
gene_loc <- gene_loc %>%
  mutate(
    START_EXT = pmax(0, START - WINDOW_SIZE),
    STOP_EXT = STOP + WINDOW_SIZE
  )

head(gene_loc)

# ==============================================================================
# 2. Load GWAS Results
# ==============================================================================

cat("\nLoading GWAS results...\n")
gwas_raw <- fread(GWAS_FILE)

cat("Total SNPs in GWAS:", format(nrow(gwas_raw), big.mark = ","), "\n")

# Rename columns for consistency
names(gwas_raw)[names(gwas_raw) == "BP"] <- "POS"

# Check p-value distribution
cat("\nGWAS p-value summary:\n")
print(summary(gwas_raw$P))

# Count significant SNPs at different thresholds
cat("\nSNPs at different thresholds:\n")
cat("  p < 1e-8:", format(sum(gwas_raw$P < 1e-8, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 1e-7:", format(sum(gwas_raw$P < 1e-7, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 1e-6:", format(sum(gwas_raw$P < 1e-6, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 1e-5:", format(sum(gwas_raw$P < 1e-5, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 1e-4:", format(sum(gwas_raw$P < 1e-4, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 1e-3:", format(sum(gwas_raw$P < 1e-3, na.rm = TRUE), big.mark = ","), "\n")
cat("  p < 0.05:", format(sum(gwas_raw$P < 0.05, na.rm = TRUE), big.mark = ","), "\n")

# ==============================================================================
# 3. Map SNPs to Genes Using 25kb Window
# ==============================================================================

cat("\n\nMapping SNPs to genes (25kb window)...\n")
cat("This may take a few minutes for large datasets...\n")

# Convert to data.table for fast overlap joins
gwas_dt <- as.data.table(gwas_raw)
gene_dt <- as.data.table(gene_loc)

# Set keys for fast overlap joins
setkey(gene_dt, CHR, START_EXT, STOP_EXT)

# Create interval columns for GWAS (point intervals)
gwas_dt[, c("START_EXT", "STOP_EXT") := .(POS, POS)]
setkey(gwas_dt, CHR, START_EXT, STOP_EXT)

# Use foverlaps for fast interval overlap
snp_gene_map <- foverlaps(
  gwas_dt,
  gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within",
  nomatch = NULL
)

cat("\nSNP-gene mappings:", format(nrow(snp_gene_map), big.mark = ","), "\n")
cat("Unique genes with mapped SNPs:", format(length(unique(snp_gene_map$GENE)), big.mark = ","), "\n")

# ==============================================================================
# 4. Aggregate to Gene-Level Statistics
# ==============================================================================

cat("\nAggregating to gene-level statistics...\n")

# Get minimum p-value per gene (best GWAS evidence for each gene)
gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    CHR = first(CHR),
    START = first(START),
    STOP = first(STOP),
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_n_snps = n(),
    gwas_best_snp = SNP[which.min(P)],
    gwas_mean_p = mean(P, na.rm = TRUE),
    gwas_median_p = median(P, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Handle any duplicate gene names
  arrange(gwas_min_p) %>%
  filter(!duplicated(GENE)) %>%
  # Sort by chromosome and position
  arrange(CHR, START)

cat("\nGenes with GWAS evidence:", format(nrow(gwas_gene), big.mark = ","), "\n")

# Genes at different significance thresholds (using min SNP p-value)
cat("\nGenes with significant GWAS hit (min SNP p-value):\n")
cat("  p < 1e-8:", sum(gwas_gene$gwas_min_p < 1e-8), "\n")
cat("  p < 1e-7:", sum(gwas_gene$gwas_min_p < 1e-7), "\n")
cat("  p < 1e-6:", sum(gwas_gene$gwas_min_p < 1e-6), "\n")
cat("  p < 1e-5:", sum(gwas_gene$gwas_min_p < 1e-5), "\n")
cat("  p < 1e-4:", sum(gwas_gene$gwas_min_p < 1e-4), "\n")
cat("  p < 1e-3:", sum(gwas_gene$gwas_min_p < 1e-3), "\n")

# Top genes
cat("\nTop 20 genes by GWAS min p-value:\n")
print(head(gwas_gene, 20))

# ==============================================================================
# 5. Save Results
# ==============================================================================

# Save gene-level results
out_file <- file.path(OUT_DIR, "AT_TN_0_5_mean_gene_results.csv")
fwrite(gwas_gene, out_file)
cat("\nGene-level results saved to:", out_file, "\n")

# Save SNP-gene mapping (optional, for reference)
out_map_file <- file.path(OUT_DIR, "AT_TN_0_5_mean_snp_gene_map.csv")
snp_gene_map_out <- snp_gene_map %>%
  select(SNP, CHR, POS, P, GENE, START, STOP)
fwrite(snp_gene_map_out, out_map_file)
cat("SNP-gene mapping saved to:", out_map_file, "\n")

# ==============================================================================
# 6. Summary Statistics
# ==============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("SUMMARY: Arabidopsis Total Nitrogen GWAS Annotation\n")
cat(strrep("=", 60), "\n")
cat("Phenotype: nitrogen_0_5cm_mean_1000\n")
cat("Window size: ", WINDOW_SIZE / 1000, "kb upstream/downstream\n", sep = "")
cat("\n")
cat("Total SNPs:", format(nrow(gwas_raw), big.mark = ","), "\n")
cat("SNPs mapped to genes:", format(nrow(snp_gene_map), big.mark = ","), "\n")
cat("Unique genes:", format(nrow(gwas_gene), big.mark = ","), "\n")
cat("\n")
cat("Output files:\n")
cat("  - ", out_file, "\n", sep = "")
cat("  - ", out_map_file, "\n", sep = "")
cat(strrep("=", 60), "\n")

# ==============================================================================
# 7. Manhattan Plot
# ==============================================================================

cat("\nGenerating Manhattan plot...\n")

# Prepare data for qqman (needs SNP, CHR, BP, P columns)
# Sample if too large for faster plotting
if (nrow(gwas_raw) > 1e6) {
  cat("  Sampling 1M SNPs for plotting (full dataset too large)...\n")
  set.seed(123)
  # Keep all significant SNPs + random sample
  sig_snps <- gwas_raw[P < 1e-4]
  nonsig_snps <- gwas_raw[P >= 1e-4]
  sample_idx <- sample(nrow(nonsig_snps), min(1e6, nrow(nonsig_snps)))
  gwas_plot <- rbind(sig_snps, nonsig_snps[sample_idx])
} else {
  gwas_plot <- gwas_raw
}

# Rename BP back to POS if needed for qqman
if (!"BP" %in% names(gwas_plot) && "POS" %in% names(gwas_plot)) {
  names(gwas_plot)[names(gwas_plot) == "POS"] <- "BP"
}

cat("  SNPs for plotting:", format(nrow(gwas_plot), big.mark = ","), "\n")

# Manhattan plot
manhattan_file <- file.path(OUT_DIR, "manhattan_AT_TN_0_5_mean.png")
png(manhattan_file, width = 2400, height = 1200, res = 300)
manhattan(
  gwas_plot,
  chr = "CHR",
  bp = "BP",
  p = "P",
  snp = "SNP",
  main = "Arabidopsis GWAS - Total Nitrogen (0-5cm Mean)",
  col = c("steelblue", "coral"),
  suggestiveline = -log10(1e-5),
  genomewideline = -log10(1e-7),
  cex = 0.6,
  cex.axis = 0.8
)
dev.off()
cat("  Manhattan plot saved:", manhattan_file, "\n")

# QQ plot
qq_file <- file.path(OUT_DIR, "qq_AT_TN_0_5_mean.png")
png(qq_file, width = 1200, height = 1200, res = 300)
qq(gwas_raw$P, main = "Arabidopsis Total Nitrogen GWAS - QQ Plot")
dev.off()
cat("  QQ plot saved:", qq_file, "\n")

# ==============================================================================
# 8. Final Summary
# ==============================================================================

cat("\n")
cat(strrep("=", 60), "\n")
cat("ALL DONE!\n")
cat(strrep("=", 60), "\n")
cat("\nSignificant SNPs:\n")
cat("  p < 1e-7 (genome-wide):", sum(gwas_raw$P < 1e-7, na.rm = TRUE), "\n")
cat("  p < 1e-5 (suggestive):", sum(gwas_raw$P < 1e-5, na.rm = TRUE), "\n")
cat("\nOutput files:\n")
cat("  - ", out_file, "\n", sep = "")
cat("  - ", out_map_file, "\n", sep = "")
cat("  - ", manhattan_file, "\n", sep = "")
cat("  - ", qq_file, "\n", sep = "")
