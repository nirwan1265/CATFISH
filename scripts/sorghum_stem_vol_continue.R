#!/usr/bin/env Rscript
# Continue CATFISH pipeline from gene results (MAGMA already completed)

library(CATFISH)

# ===== CONFIGURATION =====
out_dir    <- "sorghum_catfish_results"
out_prefix <- "sorghum_stem_vol"
n_samples  <- 357

# ===== Step 3: Combine Gene Results =====
cat("\n=== Step 3: Combining gene results from all chromosomes ===\n")

files <- list.files(
  path       = file.path(out_dir, "magma_genes"),
  pattern    = paste0(out_prefix, "_chr.*\\.multi_snp_wise\\.genes\\.out$"),
  full.names = TRUE
)
cat("Found", length(files), "gene result files\n")

gene_list <- lapply(files, function(f) {
  utils::read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, comment.char = "#")
})
genes_all_raw <- do.call(rbind, gene_list)
cat("Total gene entries:", nrow(genes_all_raw), "\n")

# Ensure P column name is correct (MAGMA multi-snp-wise uses P_SNPWISE_TOP1)
if (!"P" %in% names(genes_all_raw)) {
  if ("P_SNPWISE_TOP1" %in% names(genes_all_raw)) {
    genes_all_raw$P <- genes_all_raw$P_SNPWISE_TOP1
  } else if ("P_MULTI" %in% names(genes_all_raw)) {
    genes_all_raw$P <- genes_all_raw$P_MULTI
  }
}

cat("Column names:", paste(names(genes_all_raw), collapse = ", "), "\n")

# Deduplicate by GENE: keep smallest P
genes_all_raw <- genes_all_raw[order(genes_all_raw$GENE, genes_all_raw$P), ]
genes_all <- genes_all_raw[!duplicated(genes_all_raw$GENE), ]

# Sort by chromosome and position
if (all(c("CHR", "START") %in% names(genes_all))) {
  genes_all <- genes_all[order(genes_all$CHR, genes_all$START), ]
}

cat("Combined gene results:", nrow(genes_all), "unique genes\n")
print(head(genes_all))

# Save combined results
write.table(
  genes_all,
  file      = file.path(out_dir, paste0(out_prefix, "_genes_combined.txt")),
  sep       = "\t",
  row.names = FALSE,
  quote     = FALSE
)
cat("Saved combined gene results\n")

# ===== Step 4: Individual Pathway Tests =====
cat("\n=== Step 4: Running individual pathway tests ===\n")

# ACAT
cat("Running ACAT...\n")
acat_res <- catfish_acat_pathways(
  gene_results = genes_all,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P"
)
cat("  ACAT completed:", nrow(acat_res), "pathways\n")

# Fisher
cat("Running Fisher...\n")
fisher_res <- catfish_fisher_pathways(
  gene_results = genes_all,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P"
)
cat("  Fisher completed:", nrow(fisher_res), "pathways\n")

# Stouffer Z
cat("Running Stouffer Z...\n")
stouffer_res <- catfish_stoufferZ_pathways(
  gene_results = genes_all,
  species      = "sorghum",
  gene_col     = "GENE",
  z_col        = "ZSTAT",
  alternative  = "greater"
)
cat("  Stouffer completed:", nrow(stouffer_res), "pathways\n")

# MinP
cat("Running minP...\n")
minp_res <- catfish_minp_pathways(
  gene_results = genes_all,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P"
)
cat("  minP completed:", nrow(minp_res), "pathways\n")

# Soft TFisher
cat("Running Soft TFisher Adaptive...\n")
tfisher_res <- catfish_soft_tfisher_adaptive_pathways(
  gene_results = genes_all,
  species      = "sorghum",
  gene_col     = "GENE",
  p_col        = "P",
  tau_grid     = c(0.2, 0.1, 0.05, 0.01, 0.005, 0.001)
)
cat("  TFisher completed:", nrow(tfisher_res), "pathways\n")

# ===== Step 5: Generate Gene Correlations for MVN =====
cat("\n=== Step 5: Extracting gene correlations for MVN resampling ===\n")

out_pairs <- file.path(out_dir, paste0(out_prefix, "_gene_cor_pairs.txt"))

# Check if we have .genes.raw files for correlation extraction
raw_files <- list.files(
  path       = file.path(out_dir, "magma_genes"),
  pattern    = "\\.genes\\.raw$",
  full.names = TRUE
)

use_mvn <- FALSE  # Force global resampling for now
if (length(raw_files) > 0) {
  cat("Found", length(raw_files), "raw gene files for correlation\n")

  # Process each raw file and combine
  all_cor_pairs <- list()
  for (raw_file in raw_files) {
    cat("Processing:", basename(raw_file), "\n")
    tryCatch({
      tmp_out <- tempfile(fileext = ".txt")
      cor_df <- magma_genesraw_to_cor_pairs_banded(
        genes_raw_file = raw_file,
        out_pairs_file = tmp_out,
        gene_regex     = "^SORBI",  # Sorghum gene pattern
        keep_abs_r_ge  = 0.05,
        verbose        = FALSE
      )
      if (file.exists(tmp_out)) {
        all_cor_pairs[[basename(raw_file)]] <- data.table::fread(tmp_out)
        file.remove(tmp_out)
      }
    }, error = function(e) {
      cat("  Warning:", conditionMessage(e), "\n")
    })
  }

  if (length(all_cor_pairs) > 0) {
    combined_cor <- data.table::rbindlist(all_cor_pairs)
    data.table::fwrite(combined_cor, out_pairs, sep = "\t")
    cat("Gene correlations extracted:", nrow(combined_cor), "pairs\n")
    use_mvn <- TRUE
  }
} else {
  cat("No .genes.raw files found - using global resampling\n")
}

# ===== Step 6: Combine Results =====
cat("\n=== Step 6: Combining pathway test results ===\n")

# Combine individual results properly
results <- acat_res[, c("pathway_id", "pathway_name", "n_genes")]
results$acat_p <- acat_res$acat_p
results$fisher_p <- fisher_res$fisher_p[match(results$pathway_id, fisher_res$pathway_id)]
results$stouffer_p <- stouffer_res$stouffer_p[match(results$pathway_id, stouffer_res$pathway_id)]
results$minp_p <- minp_res$minp_p[match(results$pathway_id, minp_res$pathway_id)]
results$tfisher_p <- tfisher_res$tfisher_p_omni[match(results$pathway_id, tfisher_res$pathway_id)]

# Calculate omnibus p-value using ACAT combination of component p-values
p_cols <- c("acat_p", "fisher_p", "stouffer_p", "minp_p", "tfisher_p")
results$omni_p <- apply(results[, p_cols], 1, function(ps) {
  ps <- ps[!is.na(ps) & ps > 0 & ps < 1]
  if (length(ps) == 0) return(NA)
  # ACAT combination
  stat <- sum(tan((0.5 - ps) * pi)) / length(ps)
  0.5 - atan(stat) / pi
})

# Add FDR
results$FDR <- p.adjust(results$omni_p, method = "BH")

# Sort by omnibus p-value
results <- results[order(results$omni_p), ]

# Save final results
final_file <- file.path(out_dir, paste0(out_prefix, "_catfish_results.csv"))
write.csv(results, file = final_file, row.names = FALSE)
cat("Results saved to:", final_file, "\n")

# Show top pathways
cat("\n=== Top 20 Pathways ===\n")
print(head(results[, c("pathway_id", "pathway_name", "n_genes", "omni_p", "FDR")], 20))

# Significant pathways
sig <- results[results$FDR < 0.05, ]
cat("\n=== Significant Pathways (FDR < 0.05):", nrow(sig), "===\n")
if (nrow(sig) > 0) {
  print(sig[, c("pathway_id", "pathway_name", "n_genes", "omni_p", "FDR")])
}

cat("\n=== CATFISH Pipeline Complete ===\n")
