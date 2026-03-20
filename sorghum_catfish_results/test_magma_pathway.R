# Test MAGMA competitive pathway analysis with sorghum data

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

# Run MAGMA competitive gene-set analysis
cat("\n=== Running MAGMA competitive pathway analysis ===\n")
t1 <- Sys.time()

magma_results <- magma_geneset_competitive(
  gene_results_raw = "sorghum_catfish_results/magma_genes/sorghum_stem_vol_merged.genes.raw",
  set_annot = "sorghum_catfish_results/annot/sorghum_pathways.set.annot",
  out_prefix = "sorghum_stem_vol_magma_pathway",
  out_dir = "sorghum_catfish_results",
  pathways = pathways,
  genes_all = gene_results,
  write_tidy = TRUE
)

cat("Time:", difftime(Sys.time(), t1, units = "secs"), "seconds\n")

# Show top pathways
cat("\n=== Top 20 MAGMA pathways ===\n")
magma_results <- magma_results[order(magma_results$magma_pvalue), ]
print(head(magma_results[, c("pathway_name", "n_genes", "magma_pvalue")], 20))

# Summary
cat("\n=== Summary ===\n")
cat("Total pathways:", nrow(magma_results), "\n")
cat("Significant (p < 0.05):", sum(magma_results$magma_pvalue < 0.05, na.rm = TRUE), "\n")
cat("Significant (p < 0.01):", sum(magma_results$magma_pvalue < 0.01, na.rm = TRUE), "\n")
cat("Significant (p < 0.001):", sum(magma_results$magma_pvalue < 0.001, na.rm = TRUE), "\n")

# Save results
write.csv(magma_results, "sorghum_catfish_results/sorghum_stem_vol_MAGMA_pathways.csv", row.names = FALSE)
cat("\nMAGMA results saved to: sorghum_catfish_results/sorghum_stem_vol_MAGMA_pathways.csv\n")
