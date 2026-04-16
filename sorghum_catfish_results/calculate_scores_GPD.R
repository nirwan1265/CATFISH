# Calculate multi-layer candidate gene scores using GPD pathway results
# Combines GWAS, MAGMA, and pathway layers

library(data.table)

# ==============================================================================
# CONFIGURATION - Change tau option here
# ==============================================================================
TAU_OPTION <- "strict"  # "default" or "strict"

if (TAU_OPTION == "strict") {
  PATHWAY_FILE <- "sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD_strict_tau.csv"
  OUT_SUFFIX <- "_strict_tau"
} else {
  PATHWAY_FILE <- "sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv"
  OUT_SUFFIX <- ""
}
# ==============================================================================

# Load gene results (GWAS + MAGMA)
cat("Loading gene results...\n")
gene_results <- read.table(
  "sorghum_catfish_results/sorghum_stem_vol_genes_combined.txt",
  header = TRUE, stringsAsFactors = FALSE
)

# Load GPD pathway results (B=1M)
cat("Loading GPD pathway results (B=1M)...\n")
cat("Using:", PATHWAY_FILE, "\n")
pathway_results <- read.csv(PATHWAY_FILE, stringsAsFactors = FALSE)

# Load pathway-gene mapping
cat("Loading pathway-gene mapping...\n")
pathway_file <- "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"
pathway_raw <- read.delim(pathway_file, header = TRUE, stringsAsFactors = FALSE)
pathway_genes <- data.frame(
  pathway_id = pathway_raw[["Pathway.id"]],
  gene_id = pathway_raw[["Gene.name"]],
  stringsAsFactors = FALSE
)
pathway_genes <- unique(pathway_genes)

# Add pathway names
pathway_names <- unique(pathway_raw[, c("Pathway.id", "Pathway.name")])
names(pathway_names) <- c("pathway_id", "pathway_name")

# Calculate GWAS and MAGMA ranks
cat("Calculating ranks...\n")
gene_results$gwas_rank <- rank(gene_results$P)
gene_results$magma_rank <- rank(gene_results$P)  # Same column for this data

# Define significance thresholds
gwas_thresh <- 5e-8
magma_thresh <- 0.05 / nrow(gene_results)  # Bonferroni
pathway_thresh <- 0.05 / nrow(pathway_results)  # Bonferroni for pathways

cat("Thresholds:\n")
cat("  GWAS:", gwas_thresh, "\n")
cat("  MAGMA:", magma_thresh, "\n
")
cat("  Pathway:", pathway_thresh, "\n")

# Get top pathways (using GPD p-values)
top_pathways <- pathway_results[pathway_results$omni_p_final < pathway_thresh, ]
cat("\nTop pathways (p <", pathway_thresh, "):", nrow(top_pathways), "\n")

# For each gene, find which top pathways it belongs to
cat("\nMapping genes to top pathways...\n")
gene_pathway_info <- lapply(gene_results$GENE, function(g) {
  # Find pathways containing this gene
  gene_pws <- pathway_genes$pathway_id[pathway_genes$gene_id == g]

  # Filter to top pathways
  top_pws <- top_pathways[top_pathways$pathway_id %in% gene_pws, ]

  if (nrow(top_pws) == 0) {
    return(data.frame(
      n_top_pathways = 0,
      best_pathway_p = NA,
      pathway_rank = NA,
      pathways = NA,
      stringsAsFactors = FALSE
    ))
  }

  # Get best pathway
  best_idx <- which.min(top_pws$omni_p_final)

  data.frame(
    n_top_pathways = nrow(top_pws),
    best_pathway_p = top_pws$omni_p_final[best_idx],
    pathway_rank = rank(pathway_results$omni_p_final)[match(top_pws$pathway_id[best_idx], pathway_results$pathway_id)],
    pathways = paste(top_pws$pathway_id, collapse = "; "),
    stringsAsFactors = FALSE
  )
})
gene_pathway_df <- do.call(rbind, gene_pathway_info)

# Combine with gene results
candidates <- cbind(gene_results[, c("GENE", "P", "ZSTAT")], gene_pathway_df)
names(candidates)[2:3] <- c("magma_p", "magma_z")

# Add GWAS info (using same P column for now)
candidates$gwas_min_p <- gene_results$P
candidates$gwas_rank <- gene_results$gwas_rank
candidates$magma_rank <- gene_results$magma_rank

# Determine layer hits
candidates$hit_gwas <- candidates$gwas_min_p < gwas_thresh
candidates$hit_magma <- candidates$magma_p < magma_thresh
candidates$hit_pathway <- !is.na(candidates$n_top_pathways) & candidates$n_top_pathways > 0

# Count support layers
candidates$support_layers <- candidates$hit_gwas + candidates$hit_magma + candidates$hit_pathway

# Calculate score: -log10(gwas_rank) - log10(magma_rank) + pathway_bonus
# Higher score = better candidate
n_genes <- nrow(candidates)
candidates$score <- -log10(candidates$gwas_rank / n_genes) -
                    log10(candidates$magma_rank / n_genes)

# Add pathway bonus for genes in significant pathways
# Bonus based on -log10 of best pathway p-value (normalized)
pathway_bonus <- ifelse(
  candidates$hit_pathway,
  -log10(candidates$best_pathway_p) / 10,  # Scale down pathway contribution
  0
)
candidates$score <- candidates$score + pathway_bonus

# Sort by score
candidates <- candidates[order(-candidates$score), ]

# Show top candidates
cat("\n=== Top 20 Candidate Genes (GPD-based) ===\n")
print(head(candidates[, c("GENE", "gwas_min_p", "gwas_rank", "magma_p", "magma_rank",
                          "n_top_pathways", "best_pathway_p", "support_layers", "score")], 20))

# Summary
cat("\n=== Summary ===\n")
cat("Total genes:", nrow(candidates), "\n")
cat("Genes with 3-layer support:", sum(candidates$support_layers == 3), "\n")
cat("Genes with 2-layer support:", sum(candidates$support_layers == 2), "\n")
cat("Genes in top pathways:", sum(candidates$hit_pathway), "\n")

# Compare with old results (empirical)
cat("\n=== Comparison: Old (empirical) vs New (GPD) ===\n")
old_results <- read.csv("sorghum_catfish_results/candidate_genes_top200_sorghum_stem_vol.csv",
                        stringsAsFactors = FALSE)
cat("Old top gene:", old_results$GENE[1], "score:", old_results$score[1], "\n")
cat("New top gene:", candidates$GENE[1], "score:", candidates$score[1], "\n")

# Save results
outfile_all <- paste0("sorghum_catfish_results/candidate_genes_GPD_B1M_scored", OUT_SUFFIX, ".csv")
write.csv(candidates, outfile_all, row.names = FALSE)
cat("\nResults saved to:", outfile_all, "\n")

# Save top 200
outfile_top <- paste0("sorghum_catfish_results/candidate_genes_top200_GPD_B1M", OUT_SUFFIX, ".csv")
write.csv(head(candidates, 200), outfile_top, row.names = FALSE)
cat("Top 200 saved to:", outfile_top, "\n")
