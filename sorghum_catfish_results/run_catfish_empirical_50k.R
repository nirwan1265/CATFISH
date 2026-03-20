# Run CATFISH with EMPIRICAL mode (pure Monte Carlo) B=50,000

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

B_perm <- 50000

cat("\n========================================\n")
cat("EMPIRICAL MODE (pure Monte Carlo)\n")
cat("B =", B_perm, "permutations\n")
cat("Min possible p-value = 1/(B+1) =", format(1/(B_perm+1), scientific = TRUE), "\n")
cat("========================================\n")
cat("Start time:", format(Sys.time(), "%H:%M:%S"), "\n\n")

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
  tail_mode = "empirical",  # <-- Pure Monte Carlo, no GPD
  min_p = 1e-50
)

elapsed <- difftime(Sys.time(), t1, units = "mins")
cat("\n========================================\n")
cat("COMPLETED!\n")
cat("Total time:", round(as.numeric(elapsed), 2), "minutes\n")
cat("========================================\n")

# Show top pathways
cat("\n=== Top 20 pathways ===\n")
res <- res[order(res$omni_p_final), ]
print(head(res[, c("pathway_name", "n_genes", "omni_p_final")], 20))

# Count pathways at floor
floor_p <- 1/(B_perm + 1)
n_at_floor <- sum(res$omni_p_final <= floor_p * 1.01, na.rm = TRUE)
cat("\nPathways at permutation floor:", n_at_floor, "\n")

# Save results
outfile <- paste0("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B", B_perm, "_empirical.csv")
write.csv(res, outfile, row.names = FALSE)
cat("Results saved to:", outfile, "\n")
