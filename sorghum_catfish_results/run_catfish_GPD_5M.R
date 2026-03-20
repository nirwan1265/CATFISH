# Run CATFISH with GPD mode, B=5,000,000 (5 million)

library(devtools)
load_all(".")

# Load data
gene_results <- read.table(
  "sorghum_catfish_results/sorghum_stem_vol_genes_combined.txt",
  header = TRUE, stringsAsFactors = FALSE
)

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

cat("Loading correlation pairs...\n")
cor_pairs <- data.table::fread(
  "sorghum_catfish_results/magma_genes/magma_gene_cor_pairs_MLM_Sorghum_stem_volume.txt",
  header = TRUE, stringsAsFactors = FALSE
)
cor_pairs <- as.data.frame(cor_pairs)

B_perm <- 5000000  # 5 million
n_threads <- 12    # Use 12 of 14 cores

cat("\n============================================================\n")
cat("GPD MODE (Knijnenburg tail extrapolation) - PARALLEL\n")
cat("B =", format(B_perm, big.mark=","), "permutations\n")
cat("Threads:", n_threads, "\n")
cat("============================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

set.seed(42)
t1 <- Sys.time()

res <- catfish_omni2_pathways(
  gene_results = gene_results,
  pathways = pathways,
  gene_col = "GENE",
  p_raw_col = "P",
  z_col = "ZSTAT",
  B_perm = B_perm,
  perm_mode = "mvn",
  magma_cor_pairs = cor_pairs,
  tail_mode = "hybrid_gpd",
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  tail_min_B = 10000L,
  tail_min_tail = 50L,
  min_p = 1e-100,  # Allow even smaller p-values
  n_threads = n_threads
)

elapsed <- difftime(Sys.time(), t1, units = "hours")
cat("\n============================================================\n")
cat("COMPLETED!\n")
cat("Total time:", round(as.numeric(elapsed), 2), "hours\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================\n")

# Show top pathways
cat("\n=== Top 20 pathways ===\n")
res <- res[order(res$omni_p_final), ]
print(head(res[, c("pathway_name", "n_genes", "omni_p_final")], 20))

# Compare with B=1M results
cat("\n=== Comparison: B=1M vs B=5M ===\n")
old_res <- read.csv("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv",
                    stringsAsFactors = FALSE)
old_res <- old_res[order(old_res$omni_p_final), ]

comparison <- merge(
  data.frame(pathway_id = res$pathway_id, p_5M = res$omni_p_final, rank_5M = seq_len(nrow(res))),
  data.frame(pathway_id = old_res$pathway_id, p_1M = old_res$omni_p_final, rank_1M = seq_len(nrow(old_res))),
  by = "pathway_id"
)
comparison <- comparison[order(comparison$p_5M), ]
comparison$rank_diff <- comparison$rank_1M - comparison$rank_5M

cat("\nTop 20 pathways comparison (B=5M vs B=1M):\n")
print(head(comparison[, c("pathway_id", "p_5M", "p_1M", "rank_5M", "rank_1M", "rank_diff")], 20))

cat("\nRank correlation (Spearman):", cor(comparison$rank_5M, comparison$rank_1M, method = "spearman"), "\n")

# Save results
outfile <- "sorghum_catfish_results/sorghum_stem_vol_CATFISH_B5000000_GPD.csv"
write.csv(res, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")

# Save comparison
compfile <- "sorghum_catfish_results/comparison_B1M_vs_B5M.csv"
write.csv(comparison, compfile, row.names = FALSE)
cat("Comparison saved to:", compfile, "\n")
