#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

args_all <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_all, value = TRUE)
this_file <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "Final/final_scripts/make_supp_tables_sap_stem_volume.R"
script_dir <- dirname(normalizePath(this_file, mustWork = FALSE))
repo_dir <- normalizePath(file.path(script_dir, "..", ".."), mustWork = FALSE)

supp_dir <- file.path(repo_dir, "Final", "final_supp_tables")
results_dir <- file.path(repo_dir, "Final", "final_main_fig", "SAP_stem_volume_realdata")

gwas_file <- "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt"
gene_loc_file <- file.path(repo_dir, "inst", "extdata", "sorghum.genes.loc")
magma_file <- file.path(results_dir, "sorghum_stem_vol_genes_combined.txt")
pathway_file <- file.path(results_dir, "sorghum_stem_vol_catfish_results_MVN.csv")

WINDOW_SIZE <- 25000
GWAS_TOP_PCT <- 1
MAGMA_TOP_PCT <- 5
PATHWAY_FDR_THRESHOLD <- 0.05

dir.create(supp_dir, recursive = TRUE, showWarnings = FALSE)

write_me <- function(df, filename) {
  readr::write_csv(df, file.path(supp_dir, filename))
}

rank_seq <- function(x) {
  order_idx <- order(x, na.last = TRUE)
  out <- rep(NA_integer_, length(x))
  out[order_idx] <- seq_along(order_idx)
  out
}

pathway_name_for_gene <- function(pathway_df, gene_id) {
  hit <- pathway_df[grepl(paste0("(^|;\\s*)", gene_id, "(\\s*;|$)"), pathway_df$genes_used), , drop = FALSE]
  if (!nrow(hit)) return(NA_character_)
  hit <- hit[order(hit$omni_p_final, hit$pathway_id), , drop = FALSE]
  hit$pathway_name[[1]]
}

gene_loc <- fread(gene_loc_file, data.table = FALSE) %>%
  mutate(
    START_EXT = pmax(0, START - WINDOW_SIZE),
    STOP_EXT = STOP + WINDOW_SIZE
  )

gwas_raw <- fread(gwas_file, data.table = FALSE)
names(gwas_raw) <- sub("^chr$", "CHR", names(gwas_raw))
names(gwas_raw) <- sub("^rs$", "SNP_ID", names(gwas_raw))
names(gwas_raw) <- sub("^ps$", "POS", names(gwas_raw))
names(gwas_raw) <- sub("^p_wald$", "P", names(gwas_raw))

gwas_top50 <- gwas_raw %>%
  filter(is.finite(P), !is.na(P), P > 0) %>%
  arrange(P) %>%
  mutate(
    gwas_rank = row_number(),
    logP = -log10(P)
  ) %>%
  select(any_of(c("CHR", "SNP_ID", "POS", "P", "n_miss", "gwas_rank", "logP"))) %>%
  slice_head(n = 50)

gwas_dt <- as.data.table(gwas_raw)
gene_dt <- as.data.table(gene_loc)
setkey(gene_dt, CHR, START_EXT, STOP_EXT)
gwas_dt[, c("START_EXT", "STOP_EXT") := .(POS, POS)]
setkey(gwas_dt, CHR, START_EXT, STOP_EXT)

snp_gene_map <- foverlaps(
  gwas_dt,
  gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within",
  nomatch = NULL
)

gwas_gene <- as.data.frame(snp_gene_map) %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_n_snps = n(),
    .groups = "drop"
  ) %>%
  arrange(gwas_min_p) %>%
  filter(!duplicated(GENE))

gwas_n_top <- ceiling(nrow(gwas_gene) * GWAS_TOP_PCT / 100)
gwas_eff_threshold <- sort(gwas_gene$gwas_min_p)[gwas_n_top]
gwas_gene$gwas_rank <- rank(gwas_gene$gwas_min_p, ties.method = "min")

magma_df <- read.delim(magma_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  rename(
    magma_z = ZSTAT,
    magma_p = P
  ) %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_bonf = p.adjust(magma_p, method = "bonferroni")
  ) %>%
  arrange(magma_p)

magma_n_top <- ceiling(nrow(magma_df) * MAGMA_TOP_PCT / 100)
magma_eff_threshold <- sort(magma_df$magma_p)[magma_n_top]

magma_top50 <- magma_df %>%
  mutate(
    magma_rank = row_number(),
    logP = -log10(magma_p)
  ) %>%
  transmute(
    GENE, CHR, START, STOP, NSNPS, NPARAM, N,
    ZSTAT = magma_z, P_MULTI, P_SNPWISE_MEAN, P_SNPWISE_TOP1,
    P = magma_p, magma_rank, logP
  ) %>%
  slice_head(n = 50)

pathway_df <- read.csv(pathway_file, stringsAsFactors = FALSE)
if (!"FDR_BH" %in% names(pathway_df)) {
  if ("omni_p_final_BH" %in% names(pathway_df)) {
    pathway_df$FDR_BH <- pathway_df$omni_p_final_BH
  } else {
    pathway_df$FDR_BH <- p.adjust(pathway_df$omni_p_final, method = "BH")
  }
}

pathway_df <- pathway_df %>%
  arrange(omni_p_final, pathway_name) %>%
  mutate(pathway_rank = row_number())

top_pathways <- pathway_df %>%
  filter(FDR_BH < PATHWAY_FDR_THRESHOLD) %>%
  arrange(omni_p_final, pathway_name)

pathway_top50 <- pathway_df %>%
  transmute(
    pathway_rank,
    pathway_id,
    pathway_name,
    n_genes,
    omni_p_final,
    FDR_BH,
    dominant_component,
    agreement_score,
    genes_used
  ) %>%
  slice_head(n = 50)

pathway_genes_list <- strsplit(top_pathways$genes_used, ";")
pathway_genes <- data.frame(
  pathway_id = rep(top_pathways$pathway_id, lengths(pathway_genes_list)),
  pathway_name = rep(top_pathways$pathway_name, lengths(pathway_genes_list)),
  pathway_p = rep(top_pathways$omni_p_final, lengths(pathway_genes_list)),
  GENE = trimws(unlist(pathway_genes_list)),
  stringsAsFactors = FALSE
)

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(
    n_top_pathways = n_distinct(pathway_id),
    best_pathway_p = min(pathway_p, na.rm = TRUE),
    mean_pathway_mlog10p = mean(-log10(pathway_p), na.rm = TRUE),
    pathways = paste(pathway_id, collapse = "; "),
    .groups = "drop"
  )

gene_evidence <- magma_df %>%
  left_join(gwas_gene, by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas = !is.na(gwas_min_p) & gwas_min_p <= gwas_eff_threshold,
    hit_magma = !is.na(magma_p) & magma_p <= magma_eff_threshold,
    hit_pathway = !is.na(n_top_pathways) & n_top_pathways >= 1,
    support_layers = hit_gwas + hit_magma + hit_pathway,
    score = (hit_gwas * 1) + (hit_magma * 1) + (hit_pathway * 1) +
      0.2 * ifelse(!is.na(magma_p), -log10(magma_p), 0) +
      0.1 * ifelse(!is.na(gwas_min_p), -log10(gwas_min_p), 0) +
      0.1 * ifelse(!is.na(best_pathway_p), -log10(best_pathway_p), 0),
    magma_rank = min_rank(magma_p),
    gwas_rank = ifelse(is.na(gwas_min_p), NA_integer_, min_rank(gwas_min_p)),
    pathway_rank = ifelse(is.na(best_pathway_p), NA_integer_, min_rank(best_pathway_p)),
    evidence_class = case_when(
      hit_gwas & hit_magma & hit_pathway ~ "All 3",
      hit_gwas & hit_pathway & !hit_magma ~ "GWAS + Pathway",
      hit_magma & hit_pathway & !hit_gwas ~ "MAGMA + Pathway",
      hit_gwas & hit_magma & !hit_pathway ~ "GWAS + MAGMA",
      hit_gwas & !hit_magma & !hit_pathway ~ "GWAS only",
      hit_magma & !hit_gwas & !hit_pathway ~ "MAGMA only",
      hit_pathway & !hit_gwas & !hit_magma ~ "Pathway only",
      TRUE ~ "None"
    )
  ) %>%
  arrange(desc(score), GENE) %>%
  mutate(candidate_rank = row_number()) %>%
  transmute(
    candidate_rank,
    GENE,
    CHR,
    START,
    STOP,
    NSNPS,
    NPARAM,
    N,
    magma_z,
    P_MULTI,
    P_SNPWISE_MEAN,
    P_SNPWISE_TOP1,
    magma_p,
    magma_fdr,
    magma_bonf,
    magma_rank,
    gwas_min_p,
    gwas_n_snps,
    gwas_rank,
    n_top_pathways,
    best_pathway_p,
    mean_pathway_mlog10p,
    pathways,
    pathway_rank,
    hit_gwas,
    hit_magma,
    hit_pathway,
    support_layers,
    score,
    evidence_class
  )

component_results <- pathway_df %>%
  mutate(
    acat_rank = rank_seq(acat_p),
    fisher_rank = rank_seq(fisher_p),
    tfisher_rank = rank_seq(tfisher_p_analytic),
    minp_rank = rank_seq(minp_p_analytic),
    stouffer_rank = rank_seq(stouffer_p_analytic),
    omnibus_analytic_rank = rank_seq(omni_p_analytic),
    omnibus_final_rank = rank_seq(omni_p_final)
  ) %>%
  arrange(omni_p_final, pathway_name) %>%
  transmute(
    pathway_id,
    pathway_name,
    n_genes,
    acat_p,
    acat_rank,
    fisher_p,
    fisher_rank,
    tfisher_p_analytic,
    tfisher_rank,
    minp_p_analytic,
    minp_rank,
    stouffer_p_analytic,
    stouffer_rank,
    omni_p_analytic,
    omnibus_analytic_rank,
    omni_p_final,
    omnibus_final_rank,
    FDR_BH,
    dominant_component,
    agreement_score,
    genes_used
  )

write_me(gwas_top50, "SuppTable_7_top50_GWAS_SNPs_SAP_Stem_volume.csv")
write_me(magma_top50, "SuppTable_8_top50_MAGMA_genes_SAP_Stem_volume.csv")
write_me(pathway_top50, "SuppTable_9_top50_pathways_SAP_Stem_volume.csv")
write_me(gene_evidence, "SuppTable_10_all_candidate_genes_with_evidence_SAP_Stem_volume.csv")
write_me(component_results, "SuppTable_11_component_test_pathway_results_SAP_Stem_volume.csv")
write_me(head(gene_evidence, 200), "SuppTable_SAPStemVolume_top200_candidate_genes.csv")

message("Wrote SAP stem-volume supplementary tables to: ", supp_dir)
