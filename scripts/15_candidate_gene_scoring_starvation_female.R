#!/usr/bin/env Rscript
# =============================================================================
# Candidate Gene Scoring for Drosophila Starvation - Female
# =============================================================================
# Score = 1{GWAS} + 1{MAGMA} + 1{PATH} + 0.2*(-log10(p_MAGMA)) + 0.1*(-log10(p_GWAS))
# Thresholds: all top 1%
# =============================================================================

library(data.table)
library(dplyr)
library(tidyr)

base_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
setwd(base_dir)

message("========== STEP 1: LOAD DATA ==========")

# Load GWAS SNP results
gwas_snps <- fread("/Users/nirwantandukar/Documents/Research/data/DGRP/Starvation_stress/raw_gwas/raw_GWAS_Starvation_stress_female_DGRP.csv")
setnames(gwas_snps, "P-Value", "P.value")
setnames(gwas_snps, "Positions", "Pos")
message(sprintf("  Loaded %d SNP results", nrow(gwas_snps)))

# Load gene location
gene_loc <- fread("/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/fly.genes.loc")
message(sprintf("  Loaded %d gene locations", nrow(gene_loc)))

# Load MAGMA gene results (combine all chromosomes)
magma_files <- list.files(
  path = "/Users/nirwantandukar/Documents/Research/results/DGRP/MAGMA/Fly_magma_genes_by_chr_female",
  pattern = "\\.genes\\.out$",
  full.names = TRUE
)
magma_list <- lapply(magma_files, function(f) {
  fread(f, skip = 1)  # Skip comment line
})
gene_results <- rbindlist(magma_list)
setnames(gene_results, "P_MULTI", "P")
message(sprintf("  Loaded %d MAGMA gene results", nrow(gene_results)))

# Load CATFISH pathway results
pathway_results <- fread("catfish_omnibus_results_Fly_female_B10000/omni_acat_mvn.csv")
message(sprintf("  Loaded %d pathway results", nrow(pathway_results)))

message("\n========== STEP 2: MAP SNPS TO GENES ==========")

# Chromosome mapping: GWAS uses 2L/2R/3L/3R/X, gene_loc uses 1-5
# X=1, 2L=2, 2R=2, 3L=3, 3R=3, 4=4
gwas_snps <- gwas_snps %>%
  mutate(Chr_num = case_when(
    CHR == "X" ~ 1L,
    CHR == "2L" ~ 2L,
    CHR == "2R" ~ 2L,
    CHR == "3L" ~ 3L,
    CHR == "3R" ~ 3L,
    CHR == "4" ~ 4L,
    TRUE ~ NA_integer_
  ))

# Merge gene_loc with gene_results to get coordinates
gene_results <- gene_results %>%
  left_join(gene_loc, by = "GENE")

# For each gene, find SNPs within its boundaries
gene_gwas_info <- list()
for (i in seq_len(nrow(gene_results))) {
  g <- gene_results[i, ]
  if (is.na(g$CHR.y)) {
    gene_gwas_info[[i]] <- data.frame(GENE = g$GENE, gwas_min_p = NA, n_snps_gwas = 0)
    next
  }
  snps_in_gene <- gwas_snps %>% filter(Chr_num == g$CHR.y & Pos >= g$START.y & Pos <= g$STOP.y)
  if (nrow(snps_in_gene) > 0) {
    gene_gwas_info[[i]] <- data.frame(GENE = g$GENE, gwas_min_p = min(snps_in_gene$P.value), n_snps_gwas = nrow(snps_in_gene))
  } else {
    gene_gwas_info[[i]] <- data.frame(GENE = g$GENE, gwas_min_p = NA, n_snps_gwas = 0)
  }
}
gene_gwas_df <- rbindlist(gene_gwas_info)
message(sprintf("  Genes with GWAS SNPs: %d", sum(gene_gwas_df$n_snps_gwas > 0)))

message("\n========== STEP 3: MAP GENES TO PATHWAYS ==========")

pathway_genes_list <- list()
for (i in seq_len(nrow(pathway_results))) {
  pw <- pathway_results[i, ]
  genes <- strsplit(pw$genes_used, ";")[[1]]
  pathway_genes_list[[i]] <- data.frame(
    pathway_id = pw$pathway_id,
    gene_id = genes,
    pathway_p = pw$omni_p_final,
    stringsAsFactors = FALSE
  )
}
pathway_genes_df <- rbindlist(pathway_genes_list)

gene_pathway_info <- pathway_genes_df %>%
  group_by(gene_id) %>%
  summarise(
    n_pathways = n(),
    best_pathway_p = min(pathway_p),
    best_pathway = pathway_id[which.min(pathway_p)],
    all_pathways = paste(unique(pathway_id), collapse = ";"),
    .groups = "drop"
  )
message(sprintf("  Genes with pathway annotation: %d", nrow(gene_pathway_info)))

message("\n========== STEP 4: CALCULATE THRESHOLDS ==========")

n_genes <- nrow(gene_results)
n_snps <- nrow(gwas_snps)
n_pathways <- nrow(pathway_results)

gwas_threshold <- quantile(gwas_snps$P.value, 0.01, na.rm = TRUE)
magma_threshold <- quantile(gene_results$P, 0.01, na.rm = TRUE)
pathway_threshold_1pct <- quantile(pathway_results$omni_p_final, 0.01)
pathway_threshold_bonf <- 0.05 / n_pathways

message(sprintf("  GWAS top 1%%: p < %.4e | MAGMA top 1%%: p < %.4e", gwas_threshold, magma_threshold))
n_bonf_pathways <- sum(pathway_results$omni_p_final < pathway_threshold_bonf)
pathway_threshold <- ifelse(n_bonf_pathways > 0, pathway_threshold_bonf, pathway_threshold_1pct)
message(sprintf("  Pathways Bonferroni: %d, using %s", n_bonf_pathways, ifelse(n_bonf_pathways > 0, "Bonferroni", "top 1%")))

message("\n========== STEP 5: CALCULATE SCORES ==========")

scored_genes <- gene_results %>%
  left_join(gene_gwas_df, by = "GENE") %>%
  left_join(gene_pathway_info, by = c("GENE" = "gene_id")) %>%
  mutate(
    gwas_min_p = replace_na(gwas_min_p, 1),
    n_snps_gwas = replace_na(n_snps_gwas, 0),
    n_pathways = replace_na(n_pathways, 0),
    best_pathway_p = replace_na(best_pathway_p, 1),
    best_pathway = replace_na(best_pathway, ""),
    all_pathways = replace_na(all_pathways, "")
  ) %>%
  mutate(
    hit_gwas = as.integer(gwas_min_p < gwas_threshold),
    hit_magma = as.integer(P < magma_threshold),
    hit_pathway = as.integer(best_pathway_p < pathway_threshold)
  ) %>%
  mutate(
    magma_score = -log10(P),
    gwas_score = -log10(gwas_min_p),
    total_score = hit_gwas + hit_magma + hit_pathway + 0.2 * magma_score + 0.1 * gwas_score
  ) %>%
  arrange(P) %>%
  mutate(magma_rank = row_number()) %>%
  arrange(desc(total_score))

message("\n========== STEP 6: SUMMARY ==========")
message(sprintf("  GWAS hits: %d | MAGMA hits: %d | Pathway hits: %d",
                sum(scored_genes$hit_gwas), sum(scored_genes$hit_magma), sum(scored_genes$hit_pathway)))
message(sprintf("  2+ layers: %d | 3 layers: %d",
                sum(scored_genes$hit_gwas + scored_genes$hit_magma + scored_genes$hit_pathway >= 2),
                sum(scored_genes$hit_gwas + scored_genes$hit_magma + scored_genes$hit_pathway == 3)))

message("\n========== TOP 10 GENES ==========")
for (i in 1:min(10, nrow(scored_genes))) {
  g <- scored_genes[i, ]
  layers <- paste0(ifelse(g$hit_gwas == 1, "G", "-"),
                   ifelse(g$hit_magma == 1, "M", "-"),
                   ifelse(g$hit_pathway == 1, "P", "-"))
  message(sprintf("  %2d. %s | score=%.2f | %s | MAGMA=%.2e | GWAS=%.2e | %s",
                  i, g$GENE, g$total_score, layers, g$P, g$gwas_min_p,
                  ifelse(g$best_pathway == "", "-", g$best_pathway)))
}

# Save results
output_cols <- c("GENE", "CHR.y", "START.y", "STOP.y", "NSNPS", "ZSTAT", "P",
                 "magma_rank", "gwas_min_p", "n_snps_gwas",
                 "n_pathways", "best_pathway_p", "best_pathway", "all_pathways",
                 "hit_gwas", "hit_magma", "hit_pathway",
                 "magma_score", "gwas_score", "total_score")
output_cols <- intersect(output_cols, names(scored_genes))

out_dir <- "/Users/nirwantandukar/Documents/Research/results/DGRP/DGRP_chilipeppers/score_candidate"
fwrite(scored_genes[, ..output_cols], file.path(out_dir, "starvation_female_candidate_genes_scored.csv"))
fwrite(head(scored_genes, 200)[, ..output_cols], file.path(out_dir, "starvation_female_candidate_genes_top200.csv"))
multi_layer <- scored_genes %>% filter(hit_gwas + hit_magma + hit_pathway >= 2)
fwrite(multi_layer[, ..output_cols], file.path(out_dir, "starvation_female_genes_multi_layer.csv"))
message(sprintf("\n  Saved files to %s", out_dir))
