# Run CATFISH with GPD mode, B=1,000,000 for Arabidopsis
# Parallel processing enabled

library(devtools)
load_all(".")
library(dplyr)

# ============================================================
# 1. Load and prepare gene results from MAGMA
# ============================================================
cat("Loading MAGMA gene results...\n")

files <- list.files(
  path = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/MAGMA/AT_cold_by_chr",
  pattern = "^AT_cold_chr.*\\.multi_snp_wise.genes\\.out$",
  full.names = TRUE
)

gene_list <- lapply(files, function(f) {
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})

genes_all_raw <- do.call(rbind, gene_list)
colnames(genes_all_raw)[9] <- "P"

# Keep best P per gene
o <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

cat("Total genes:", nrow(genes_all), "\n")

# ============================================================
# 2. Adjust p-values based on gene length and NSNPS
# ============================================================
cat("Loading gene lengths and adjusting p-values...\n")

gene_len <- read.delim("inst/extdata/Arabidopsis_gene_lengths.tsv")
gene_len$gene_id <- sub("^gene:", "", gene_len$gene_id)

adj_out <- catfish_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = gene_len,
  gene_col     = "GENE",
  nsnp_col     = "NSNPS",
  p_col        = "P",
  z_col        = "ZSTAT",
  len_gene_col = "gene_id",
  len_col      = "length"
)

genes_adj <- adj_out[, c(1, 2, 3)]
colnames(genes_adj) <- c("GENE", "ZSTAT", "P")
cat("Adjusted genes:", nrow(genes_adj), "\n")

# ============================================================
# 3. Load pathways
# ============================================================
cat("Loading Arabidopsis pathways...\n")

at_pw <- read.delim("inst/extdata/pathway/aracyc_pathways.20230103")
at_pw <- at_pw %>%
  transmute(
    pathway_id   = Pathway.id,
    pathway_name = Pathway.name,
    gene_id      = Gene.id
  ) %>%
  mutate(
    gene_id = sub("^gene:", "", gene_id)
  ) %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  distinct()

cat("Pathways loaded:", length(unique(at_pw$pathway_id)), "\n")

# ============================================================
# 4. Load correlation pairs
# ============================================================
cat("Loading correlation pairs...\n")

cor_pairs <- data.table::fread(
  "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM_arabidopsis.txt",
  header = TRUE, stringsAsFactors = FALSE
)
cor_pairs <- as.data.frame(cor_pairs)

# ============================================================
# 5. Run CATFISH with GPD mode
# ============================================================
B_perm <- 1000000  # 1 million
n_threads <- 12    # Use 12 of 14 cores

cat("\n============================================================\n")
cat("GPD MODE (Knijnenburg tail extrapolation) - PARALLEL\n")
cat("B =", format(B_perm, big.mark = ","), "permutations\n")
cat("Threads:", n_threads, "\n")
cat("============================================================\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

set.seed(42)
t1 <- Sys.time()

res <- catfish_omni2_pathways(
  gene_results = genes_adj,
  pathways = at_pw,
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
  min_p = 1e-50,
  n_threads = n_threads,
  output = TRUE,
  out_dir = "arabidopsis_catfish_results/catfish_GPD_B1M"
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

# Save results
outfile <- "arabidopsis_catfish_results/arabidopsis_cold_CATFISH_B1000000_GPD.csv"
write.csv(res, outfile, row.names = FALSE)
cat("\nResults saved to:", outfile, "\n")

# Also save genes_adj for scoring later
write.table(genes_adj, "arabidopsis_catfish_results/arabidopsis_cold_genes_adjusted.txt",
            row.names = FALSE, quote = FALSE)
cat("Gene results saved to: arabidopsis_catfish_results/arabidopsis_cold_genes_adjusted.txt\n")
