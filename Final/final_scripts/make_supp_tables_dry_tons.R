suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

base_dir <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre"
out_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_supp_tables"

catfish_file <- file.path(
  base_dir,
  "CATFISH_permutation_B10000_mvn_GPD_paper_tau_false",
  "Dry_tons_per_acre_CATFISH_ACAT_mvn_B10000_GPD_paper_tau_false.csv"
)

candidate_file <- file.path(
  base_dir,
  "candidate_gene_scoring_B10000_GPD_paper_tau_false",
  "candidate_genes_all_Dry_tons_per_acre_B10000_GPD_paper_tau_false.csv"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

catfish_df <- read_csv(catfish_file, show_col_types = FALSE)
candidate_df <- read_csv(candidate_file, show_col_types = FALSE)

bonf_thresh <- 0.05 / nrow(catfish_df)

rank_seq <- function(x) {
  order_idx <- order(x, na.last = TRUE)
  out <- rep(NA_integer_, length(x))
  out[order_idx] <- seq_along(order_idx)
  out
}

fig6_top_pathways <- catfish_df %>%
  arrange(omni_p_final, pathway_name) %>%
  mutate(pathway_rank = row_number()) %>%
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

fig6_candidate_top50 <- candidate_df %>%
  arrange(desc(score), GENE) %>%
  mutate(candidate_rank = row_number()) %>%
  select(candidate_rank, everything()) %>%
  slice_head(n = 50)

fig6_candidate_all <- candidate_df %>%
  arrange(desc(score), GENE) %>%
  mutate(candidate_rank = row_number()) %>%
  select(candidate_rank, everything())

component_cols <- c(
  ACAT = "acat_p",
  Fisher = "fisher_p",
  TFisher = "tfisher_p_analytic",
  minP = "minp_p_analytic",
  Stouffer = "stouffer_p_analytic",
  Omnibus_analytic = "omni_p_analytic",
  Omnibus_final = "omni_p_final"
)

fig7_method_summary <- bind_rows(lapply(names(component_cols), function(method) {
  col <- component_cols[[method]]
  tibble(
    method = method,
    p_column = col,
    n_pathways_tested = sum(is.finite(catfish_df[[col]])),
    bonferroni_threshold = bonf_thresh,
    n_bonferroni_significant = sum(catfish_df[[col]] < bonf_thresh, na.rm = TRUE),
    min_p = min(catfish_df[[col]], na.rm = TRUE),
    median_p = median(catfish_df[[col]], na.rm = TRUE),
    max_p = max(catfish_df[[col]], na.rm = TRUE)
  )
}))

fig7_component_results <- catfish_df %>%
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

candidate_with_class <- fig6_candidate_all %>%
  mutate(
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
  )

fig7_candidate_summary <- candidate_with_class %>%
  count(evidence_class, name = "n_genes") %>%
  arrange(desc(n_genes), evidence_class)

readr::write_csv(
  fig6_top_pathways,
  file.path(out_dir, "SuppTable_Fig6_top50_pathways_Dry_tons_per_acre.csv")
)

readr::write_csv(
  fig6_candidate_top50,
  file.path(out_dir, "SuppTable_Fig6_top50_candidate_genes_Dry_tons_per_acre.csv")
)

readr::write_csv(
  fig6_candidate_all,
  file.path(out_dir, "SuppTable_Fig6_all_candidate_genes_Dry_tons_per_acre.csv")
)

readr::write_csv(
  fig7_method_summary,
  file.path(out_dir, "SuppTable_Fig7_method_significance_summary_Dry_tons_per_acre.csv")
)

readr::write_csv(
  fig7_component_results,
  file.path(out_dir, "SuppTable_Fig7_component_test_pathway_results_Dry_tons_per_acre.csv")
)

readr::write_csv(
  fig7_candidate_summary,
  file.path(out_dir, "SuppTable_Fig7_candidate_evidence_summary_Dry_tons_per_acre.csv")
)

readr::write_csv(
  candidate_with_class,
  file.path(out_dir, "SuppTable_Fig7_all_candidate_genes_with_evidence_Dry_tons_per_acre.csv")
)

readLines(textConnection(c(
  "Supplementary tables exported from final Dry tons per acre figure inputs.",
  "",
  "Files:",
  "- SuppTable_Fig6_top50_pathways_Dry_tons_per_acre.csv",
  "- SuppTable_Fig6_top50_candidate_genes_Dry_tons_per_acre.csv",
  "- SuppTable_Fig6_all_candidate_genes_Dry_tons_per_acre.csv",
  "- SuppTable_Fig7_method_significance_summary_Dry_tons_per_acre.csv",
  "- SuppTable_Fig7_component_test_pathway_results_Dry_tons_per_acre.csv",
  "- SuppTable_Fig7_candidate_evidence_summary_Dry_tons_per_acre.csv",
  "- SuppTable_Fig7_all_candidate_genes_with_evidence_Dry_tons_per_acre.csv"
))) %>%
  writeLines(con = file.path(out_dir, "README.txt"))

cat("Supplementary tables written to:\n", out_dir, "\n", sep = "")
