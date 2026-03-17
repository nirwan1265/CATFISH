# Test GPD tail extrapolation with sorghum stem volume data
# This script compares empirical vs hybrid_gpd tail modes

library(devtools)
load_all(".")

# Check if evir is installed
if (!requireNamespace("evir", quietly = TRUE)) {
  message("Installing evir package for GPD tail fitting...")
  install.packages("evir")
}

# Load the sorghum gene results
gene_results <- read.table(
  "sorghum_catfish_results/sorghum_stem_vol_genes_combined.txt",
  header = TRUE, stringsAsFactors = FALSE
)

# Load gene correlations (full file with ~2.7M pairs)
cat("Loading correlation pairs (this may take a moment)...\n")
cor_pairs <- data.table::fread(
  "sorghum_catfish_results/magma_genes/magma_gene_cor_pairs_MLM_Sorghum_stem_volume.txt",
  header = TRUE, stringsAsFactors = FALSE
)
cor_pairs <- as.data.frame(cor_pairs)

cat("Gene results:", nrow(gene_results), "genes\n")
cat("Correlation pairs:", nrow(cor_pairs), "pairs\n")

# Run with EMPIRICAL mode (original behavior)
cat("\n=== Running with tail_mode = 'empirical' (original behavior) ===\n")
set.seed(42)
res_empirical <- catfish_omni2_pathways(
  gene_results = gene_results,
  species = "sorghum",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = 10000,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "empirical",
  min_p = 1e-50
)

# Run with HYBRID_GPD mode (new behavior)
cat("\n=== Running with tail_mode = 'hybrid_gpd' (GPD extrapolation) ===\n")
set.seed(42)
res_gpd <- catfish_omni2_pathways(
  gene_results = gene_results,
  species = "sorghum",
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = 10000,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "hybrid_gpd",
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  tail_min_B = 10000L,
  tail_min_tail = 50L,
  min_p = 1e-50
)

# Compare results
cat("\n=== Comparison of top 20 pathways ===\n")
comparison <- merge(
  res_empirical[, c("pathway_id", "pathway_name", "omni_p_final")],
  res_gpd[, c("pathway_id", "omni_p_final")],
  by = "pathway_id",
  suffixes = c("_emp", "_gpd")
)

# Order by empirical p-value
comparison <- comparison[order(comparison$omni_p_final_emp), ]
comparison$improved <- comparison$omni_p_final_gpd < comparison$omni_p_final_emp

cat("\nTop 20 pathways:\n")
print(head(comparison[, c("pathway_name", "omni_p_final_emp", "omni_p_final_gpd", "improved")], 20))

# Check how many pathways hit the floor
n_floor_emp <- sum(comparison$omni_p_final_emp <= 1e-4, na.rm = TRUE)
n_floor_gpd <- sum(comparison$omni_p_final_gpd <= 1e-4, na.rm = TRUE)

cat("\n=== Summary ===\n")
cat("Pathways at permutation floor (p <= 1e-4):\n")
cat("  Empirical:", n_floor_emp, "\n")
cat("  GPD:", n_floor_gpd, "\n")

# Show pathways that now have smaller p-values
improved <- comparison[comparison$improved & !is.na(comparison$improved), ]
cat("\nPathways improved by GPD:", nrow(improved), "\n")

if (nrow(improved) > 0) {
  cat("\nImproved pathways (showing first 10):\n")
  print(head(improved[, c("pathway_name", "omni_p_final_emp", "omni_p_final_gpd")], 10))
}

# Save results
write.csv(res_gpd, "sorghum_catfish_results/sorghum_stem_vol_catfish_results_GPD.csv", row.names = FALSE)
cat("\nGPD results saved to: sorghum_catfish_results/sorghum_stem_vol_catfish_results_GPD.csv\n")
