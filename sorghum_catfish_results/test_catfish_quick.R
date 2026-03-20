# Quick CATFISH pathway analysis with B=100 for comparison with MAGMA

library(devtools)
load_all(".")

# Load the sorghum gene results
gene_results <- read.table(
  "sorghum_catfish_results/sorghum_stem_vol_genes_combined.txt",
  header = TRUE, stringsAsFactors = FALSE
)

# Load the SORBI-format pathway file
cat("Loading SORBI-format pathway file...\n")
pathway_file <- "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"
pathway_raw <- read.delim(pathway_file, header = TRUE, stringsAsFactors = FALSE)
pathways <- data.frame(
  pathway_id = pathway_raw[["Pathway.id"]],
  pathway_name = pathway_raw[["Pathway.name"]],
  gene_id = pathway_raw[["Gene.name"]],
  stringsAsFactors = FALSE
)
pathways <- unique(pathways)
cat("Loaded", nrow(pathways), "pathway-gene pairs\n")

# Load gene correlations
cat("Loading correlation pairs...\n")
cor_pairs <- data.table::fread(
  "sorghum_catfish_results/magma_genes/magma_gene_cor_pairs_MLM_Sorghum_stem_volume.txt",
  header = TRUE, stringsAsFactors = FALSE
)
cor_pairs <- as.data.frame(cor_pairs)

# Load MAGMA pathway results for comparison
cat("Loading MAGMA pathway results...\n")
magma_results <- read.csv("sorghum_catfish_results/sorghum_stem_vol_MAGMA_pathways.csv",
                          stringsAsFactors = FALSE)

# Run CATFISH with B=100 (quick test)
cat("\n=== Running CATFISH with B=100 ===\n")
set.seed(42)
t1 <- Sys.time()
catfish_results <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = 100,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  magma_out = magma_results[, c("pathway_id", "magma_pvalue")],  # Include MAGMA
  include_magma_in_omni = TRUE,
  tail_mode = "empirical",
  min_p = 1e-50
)
cat("Time:", difftime(Sys.time(), t1, units = "secs"), "seconds\n")

# Show top pathways
cat("\n=== Top 20 CATFISH pathways ===\n")
catfish_results <- catfish_results[order(catfish_results$omni_p_final), ]
print(head(catfish_results[, c("pathway_name", "n_genes", "omni_p_final", "magma_pvalue")], 20))

# Compare with MAGMA
cat("\n=== Comparison: CATFISH vs MAGMA ===\n")
comparison <- merge(
  catfish_results[, c("pathway_id", "pathway_name", "omni_p_final")],
  magma_results[, c("pathway_id", "magma_pvalue")],
  by = "pathway_id"
)
comparison$catfish_rank <- rank(comparison$omni_p_final)
comparison$magma_rank <- rank(comparison$magma_pvalue)
comparison <- comparison[order(comparison$omni_p_final), ]

cat("\nCorrelation between CATFISH and MAGMA p-values:\n")
cat("  Spearman:", cor(comparison$omni_p_final, comparison$magma_pvalue, method = "spearman"), "\n")
cat("  Pearson (log10):", cor(log10(comparison$omni_p_final), log10(comparison$magma_pvalue)), "\n")

# Save results
write.csv(catfish_results, "sorghum_catfish_results/sorghum_stem_vol_CATFISH_B100.csv", row.names = FALSE)
cat("\nCATFISH results saved to: sorghum_catfish_results/sorghum_stem_vol_CATFISH_B100.csv\n")
