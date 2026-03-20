# Quick test of GPD tail extrapolation with sorghum stem volume data
# Uses fewer permutations for faster testing

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

# Load the SORBI-format pathway file (matches SORBI_3xxx gene IDs in gene_results)
cat("Loading SORBI-format pathway file...\n")
pathway_file <- "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"
pathway_raw <- read.delim(pathway_file, header = TRUE, stringsAsFactors = FALSE)
# Convert to required format
pathways <- data.frame(
  pathway_id = pathway_raw[["Pathway.id"]],
  pathway_name = pathway_raw[["Pathway.name"]],
  gene_id = pathway_raw[["Gene.name"]],  # SORBI format
  stringsAsFactors = FALSE
)
pathways <- unique(pathways)
cat("Loaded", nrow(pathways), "pathway-gene pairs\n")

# Load gene correlations (full file with ~2.7M pairs)
cat("Loading correlation pairs (this may take a moment)...\n")
cor_pairs <- data.table::fread(
  "sorghum_catfish_results/magma_genes/magma_gene_cor_pairs_MLM_Sorghum_stem_volume.txt",
  header = TRUE, stringsAsFactors = FALSE
)
cor_pairs <- as.data.frame(cor_pairs)

cat("Gene results:", nrow(gene_results), "genes\n")
cat("Correlation pairs:", nrow(cor_pairs), "pairs\n")

# Quick test with fewer permutations
B_test <- 1000

# Run with EMPIRICAL mode (original behavior)
cat("\n=== Running with tail_mode = 'empirical' (B=", B_test, ") ===\n")
set.seed(42)
t1 <- Sys.time()
res_empirical <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,  # Use custom SORBI-format pathways
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = B_test,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "empirical",
  min_p = 1e-50
)
cat("Time:", difftime(Sys.time(), t1, units = "secs"), "seconds\n")

# Run with HYBRID_GPD mode (new behavior)
cat("\n=== Running with tail_mode = 'hybrid_gpd' (B=", B_test, ") ===\n")
set.seed(42)
t1 <- Sys.time()
res_gpd <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,  # Use custom SORBI-format pathways
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = B_test,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "hybrid_gpd",
  tail_switch_exceed = 10L,
  tail_gpd_k = 100L,  # smaller for quick test
  tail_min_B = 500L,  # lower threshold for quick test
  tail_min_tail = 25L,
  min_p = 1e-50
)
cat("Time:", difftime(Sys.time(), t1, units = "secs"), "seconds\n")

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
comparison$ratio <- comparison$omni_p_final_gpd / comparison$omni_p_final_emp

cat("\nTop 20 pathways (sorted by empirical p-value):\n")
print(head(comparison[, c("pathway_name", "omni_p_final_emp", "omni_p_final_gpd", "ratio")], 20))

# Check how many pathways hit the floor
floor_val <- 1 / (B_test + 1)
n_floor_emp <- sum(comparison$omni_p_final_emp <= floor_val * 1.1, na.rm = TRUE)
n_floor_gpd <- sum(comparison$omni_p_final_gpd <= floor_val * 1.1, na.rm = TRUE)

cat("\n=== Summary ===\n")
cat("Permutation floor (1/B+1):", format(floor_val, scientific = TRUE), "\n")
cat("Pathways at/near floor:\n")
cat("  Empirical:", n_floor_emp, "\n")
cat("  GPD:", n_floor_gpd, "\n")

# Show pathways where GPD gave smaller p-values
improved <- comparison[comparison$ratio < 1 & !is.na(comparison$ratio), ]
cat("\nPathways improved by GPD:", nrow(improved), "\n")

cat("\nTest completed successfully!\n")
