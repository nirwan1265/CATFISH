# Run CATFISH with progress tracking
# You can see progress and cancel if needed

library(devtools)
load_all(".")  # Reload to pick up any changes

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

cat("Gene results:", nrow(gene_results), "genes\n")
cat("Correlation pairs:", nrow(cor_pairs), "pairs\n")

# Set B value - change this to adjust permutations
B_perm <- 10000  # 10^4 permutations

cat("\n========================================\n")
cat("Starting CATFISH with B =", B_perm, "permutations\n")
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
  tail_mode = "hybrid_gpd",
  tail_switch_exceed = 10L,
  tail_gpd_k = 250L,
  tail_min_B = 10000L,
  tail_min_tail = 50L,
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

# Save results
outfile <- paste0("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B", B_perm, "_GPD.csv")
write.csv(res, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")
