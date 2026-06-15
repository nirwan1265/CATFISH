suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
  library(tidyr)
})

repo_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
final_dir <- file.path(repo_dir, "Final")
out_dir <- file.path(final_dir, "final_supp_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

write_me <- function(df, filename) {
  readr::write_csv(df, file.path(out_dir, filename))
}

rank_seq <- function(x) {
  order_idx <- order(x, na.last = TRUE)
  out <- rep(NA_integer_, length(x))
  out[order_idx] <- seq_along(order_idx)
  out
}

pick_component_col <- function(df, candidates) {
  for (col in candidates) {
    if (!col %in% names(df)) next
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (any(is.finite(vals) & !is.na(vals))) return(col)
  }
  stop("Could not find a populated component column among: ", paste(candidates, collapse = ", "), call. = FALSE)
}

# ------------------------------------------------------------------------------
# Fig1: archetype rank curves
# ------------------------------------------------------------------------------
set.seed(1)

rank_curve <- function(p, label, group = "in_set") {
  p <- p[is.finite(p) & !is.na(p)]
  p <- pmin(pmax(p, .Machine$double.xmin), 1)
  tibble(
    archetype = label,
    group = group,
    rank = seq_along(p) / length(p),
    p = sort(p),
    mlog10p = -log10(sort(p))
  )
}

make_exchangeable_sigma <- function(m, rho) {
  matrix(rho, nrow = m, ncol = m) + diag(1 - rho, m)
}

simulate_mvn_p <- function(mu, sigma, one_sided = TRUE) {
  m <- length(mu)
  z <- as.numeric(mu + t(chol(sigma)) %*% rnorm(m))
  p <- if (one_sided) pnorm(z, lower.tail = FALSE) else 2 * pnorm(-abs(z))
  pmin(pmax(p, .Machine$double.xmin), 1)
}

m <- 220
rho <- 0.20
sigma <- make_exchangeable_sigma(m = m, rho = rho)

mu_I <- rep(0, m)
mu_I[sample.int(m, 6)] <- 3.2
p_I <- simulate_mvn_p(mu_I, sigma)

mu_II <- rep(0, m)
mu_II[sample.int(m, round(0.45 * m))] <- 1.35
p_II <- simulate_mvn_p(mu_II, sigma)

mu_III <- rep(0, m)
mu_III[sample.int(m, round(0.85 * m))] <- 0.45
p_III <- simulate_mvn_p(mu_III, sigma)

mu_IV <- rep(0, m)
idx_driver <- sample.int(m, 4)
idx_support <- sample(setdiff(seq_len(m), idx_driver), 40)
mu_IV[idx_driver] <- 3.0
mu_IV[idx_support] <- 1.1
p_IV <- simulate_mvn_p(mu_IV, sigma)

mu_V <- rep(0, m)
mu_V[sample.int(m, 1)] <- 4.2
p_V <- simulate_mvn_p(mu_V, sigma)

fig1_curves <- bind_rows(
  rank_curve(p_I, "Archetype I - SDA"),
  rank_curve(p_II, "Archetype II - CME"),
  rank_curve(p_III, "Archetype III - DPS"),
  rank_curve(p_IV, "Archetype IV - HDS"),
  rank_curve(p_V, "Archetype V - SGP")
)

fig1_params <- tibble(
  parameter = c("seed", "m", "rho", "archetype_I_drivers", "archetype_II_fraction",
                "archetype_III_fraction", "archetype_IV_drivers", "archetype_IV_support",
                "archetype_V_drivers"),
  value = c("1", "220", "0.20", "6", "0.45", "0.85", "4", "40", "1")
)

write_me(fig1_curves, "SuppTable_Fig1_archetype_rank_curves.csv")
write_me(fig1_params, "SuppTable_Fig1_simulation_parameters.csv")

# ------------------------------------------------------------------------------
# Fig3: existing simulation tables
# ------------------------------------------------------------------------------
fig3_table_dir <- file.path(final_dir, "final_sup", "tables")
fig3_compact <- file.path(fig3_table_dir, "TableS_Archetype_Verification_ByG_Compact.csv")
fig3_full <- file.path(fig3_table_dir, "TableS_Archetype_Verification_ByG_Full.csv")

if (file.exists(fig3_compact)) {
  file.copy(fig3_compact, file.path(out_dir, "SuppTable_Fig3_archetype_verification_compact.csv"), overwrite = TRUE)
}
if (file.exists(fig3_full)) {
  file.copy(fig3_full, file.path(out_dir, "SuppTable_Fig3_archetype_verification_full.csv"), overwrite = TRUE)
}

block_a_csv <- file.path(repo_dir, "simulation_results", "block_a_archetype_summary_by_m.csv")
if (file.exists(block_a_csv)) {
  file.copy(block_a_csv, file.path(out_dir, "SuppTable_Fig3_source_block_a_by_m.csv"), overwrite = TRUE)
}

# ------------------------------------------------------------------------------
# Fig4: null calibration + missing-correlation calibration
# ------------------------------------------------------------------------------
sim_dir <- file.path(repo_dir, "simulation_results")
block_b <- readRDS(file.path(sim_dir, "block_b.rds"))
block_c <- readRDS(file.path(sim_dir, "block_c.rds"))

write_me(block_b, "SuppTable_Fig4_null_calibration_blockB.csv")
write_me(block_c, "SuppTable_Fig4_missing_correlation_blockC.csv")

# ------------------------------------------------------------------------------
# Fig5: component-breakage + adaptive / leave-one-out
# ------------------------------------------------------------------------------
block_d <- readRDS(file.path(sim_dir, "block_d.rds"))
block_e_adaptive <- readRDS(file.path(sim_dir, "block_e_adaptive.rds"))
block_e_leave1out <- readRDS(file.path(sim_dir, "block_e_leave1out.rds"))

write_me(block_d, "SuppTable_Fig5_component_breakage_blockD.csv")
write_me(block_e_adaptive, "SuppTable_Fig5_adaptive_blockE.csv")
write_me(block_e_leave1out, "SuppTable_Fig5_leave_one_out_blockE.csv")

adaptive_summary <- block_e_adaptive %>%
  mutate(
    source = "Adaptive",
    method_label = case_when(
      method == "omnibus_analytical" ~ "Analytic",
      method == "omnibus_mvn" ~ "MVN",
      method == "omnibus_adaptive" ~ "Adaptive+Analytic",
      method == "omnibus_adaptive_mvn" ~ "Adaptive+MVN",
      TRUE ~ method
    )
  ) %>%
  group_by(source, method_label) %>%
  summarise(
    lambda = mean(lambda, na.rm = TRUE),
    lambda_se = sd(lambda, na.rm = TRUE) / sqrt(n()),
    type1_05 = mean(type1_05, na.rm = TRUE),
    type1_05_se = sd(type1_05, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

leave1out_summary <- block_e_leave1out %>%
  mutate(
    source = "Leave-one-out",
    method_label = case_when(
      method == "omnibus_minus_acat" ~ "-ACAT",
      method == "omnibus_minus_fisher" ~ "-Fisher",
      method == "omnibus_minus_tfisher" ~ "-TFisher",
      method == "omnibus_minus_minp" ~ "-minP",
      method == "omnibus_minus_stouffer" ~ "-Stouffer",
      method == "omnibus_all" ~ "All",
      TRUE ~ method
    )
  ) %>%
  group_by(source, method_label) %>%
  summarise(
    lambda = mean(lambda, na.rm = TRUE),
    lambda_se = sd(lambda, na.rm = TRUE) / sqrt(n()),
    type1_05 = mean(type1_05, na.rm = TRUE),
    type1_05_se = sd(type1_05, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

write_me(bind_rows(adaptive_summary, leave1out_summary), "SuppTable_Fig5_adaptive_leave1out_summary.csv")

# ------------------------------------------------------------------------------
# Fig6: Dry tons GWAS + MAGMA + CATFISH
# ------------------------------------------------------------------------------
gwas_dir <- "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/BAP/Dry_tons_per_acre"
catfish_dir <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre/CATFISH_permutation_B1000000_mvn_GPD_strict_tau"
magma_file <- file.path(catfish_dir, "Dry_tons_per_acre_combined_genes.tsv")
pathway_file <- file.path(catfish_dir, "Dry_tons_per_acre_CATFISH_ACAT_mvn_B1000000_GPD_strict_tau.csv")
candidate_file <- file.path(
  "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre/candidate_gene_scoring_B1000000_GPD_strict_tau",
  "candidate_genes_all_Dry_tons_per_acre_B1000000_GPD_strict_tau.csv"
)

gwas_files <- list.files(gwas_dir, pattern = "\\.assoc\\.txt$", full.names = TRUE)
gwas_top <- bind_rows(lapply(gwas_files, function(f) {
  x <- fread(f, header = TRUE, data.table = FALSE)
  names(x) <- sub("^chr$", "CHR", names(x))
  names(x) <- sub("^rs$", "SNP_ID", names(x))
  names(x) <- sub("^ps$", "POS", names(x))
  names(x) <- sub("^p_wald$", "P", names(x))
  req <- intersect(c("CHR", "SNP_ID", "POS", "P", "n_miss"), names(x))
  x <- x[, req, drop = FALSE]
  x <- x[is.finite(x$P) & x$P > 0, , drop = FALSE]
  x <- x[order(x$P), , drop = FALSE]
  head(x, 50)
})) %>%
  arrange(P) %>%
  mutate(gwas_rank = row_number(), logP = -log10(P)) %>%
  slice_head(n = 50)

write_me(gwas_top, "SuppTable_Fig6_top50_GWAS_SNPs_Dry_tons_per_acre.csv")

magma_top <- read.delim(magma_file, stringsAsFactors = FALSE, check.names = FALSE) %>%
  arrange(P) %>%
  mutate(magma_rank = row_number(), logP = -log10(P)) %>%
  slice_head(n = 50)

write_me(magma_top, "SuppTable_Fig6_top50_MAGMA_genes_Dry_tons_per_acre.csv")

pathway_results <- fread(pathway_file, data.table = FALSE)
acat_col <- pick_component_col(pathway_results, c("acat_p_mvn_cal", "acat_p"))
fisher_col <- pick_component_col(pathway_results, c("fisher_p_mvn_cal", "fisher_p"))
stouffer_col <- pick_component_col(pathway_results, c("stouffer_p_mvn_cal", "stouffer_p_analytic"))
minp_col <- pick_component_col(pathway_results, c("minp_p_mvn_cal", "minp_p_analytic"))
tfisher_col <- pick_component_col(pathway_results, c("tfisher_p_mvn_cal", "tfisher_p_analytic"))
omni_col <- pick_component_col(pathway_results, c("omni_p_final", "omni_p_mvn", "omni_p_analytic"))

fig6_heatmap <- pathway_results %>%
  arrange(.data[[omni_col]], pathway_name) %>%
  transmute(
    pathway_rank = row_number(),
    pathway_id,
    pathway_name,
    n_genes,
    ACAT = .data[[acat_col]],
    Fisher = .data[[fisher_col]],
    Stouffer = .data[[stouffer_col]],
    minP = .data[[minp_col]],
    TFisher = .data[[tfisher_col]],
    Omnibus = .data[[omni_col]]
  ) %>%
  slice_head(n = 18)

write_me(fig6_heatmap, "SuppTable_Fig6_top18_pathway_heatmap_Dry_tons_per_acre.csv")

fig6_top50_pathways <- pathway_results %>%
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

write_me(fig6_top50_pathways, "SuppTable_Fig6_top50_pathways_Dry_tons_per_acre.csv")

# ------------------------------------------------------------------------------
# Fig7: CATFISH method comparison + integrated candidates
# ------------------------------------------------------------------------------
catfish_df <- pathway_results
candidate_df <- read_csv(candidate_file, show_col_types = FALSE)
bonf_thresh <- 0.05 / nrow(catfish_df)

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

write_me(fig7_method_summary, "SuppTable_Fig7_method_significance_summary_Dry_tons_per_acre.csv")

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

write_me(fig7_component_results, "SuppTable_Fig7_component_test_pathway_results_Dry_tons_per_acre.csv")

candidate_with_class <- candidate_df %>%
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
  ) %>%
  arrange(desc(score), GENE) %>%
  mutate(candidate_rank = row_number()) %>%
  select(candidate_rank, everything())

write_me(candidate_with_class, "SuppTable_Fig7_all_candidate_genes_with_evidence_Dry_tons_per_acre.csv")
write_me(head(candidate_with_class, 50), "SuppTable_Fig7_top50_candidate_genes_Dry_tons_per_acre.csv")

top10_gwas_genes <- candidate_with_class %>%
  filter(!is.na(gwas_rank)) %>%
  arrange(gwas_rank, GENE) %>%
  transmute(
    GENE,
    score,
    gwas_rank,
    gwas_min_p,
    magma_rank,
    magma_p,
    pathway_rank,
    best_pathway_name,
    n_top_pathways,
    hit_gwas,
    hit_magma,
    hit_pathway
  ) %>%
  slice_head(n = 10)

top10_magma_genes <- candidate_with_class %>%
  filter(!is.na(magma_rank)) %>%
  arrange(magma_rank, GENE) %>%
  transmute(
    GENE,
    score,
    magma_rank,
    magma_p,
    gwas_rank,
    gwas_min_p,
    pathway_rank,
    best_pathway_name,
    n_top_pathways,
    hit_gwas,
    hit_magma,
    hit_pathway
  ) %>%
  slice_head(n = 10)

top10_all3_genes <- candidate_with_class %>%
  filter(hit_gwas, hit_magma, hit_pathway) %>%
  arrange(desc(score), GENE) %>%
  transmute(
    GENE,
    score,
    gwas_rank,
    gwas_min_p,
    magma_rank,
    magma_p,
    pathway_rank,
    best_pathway_name,
    best_pathway_p,
    n_top_pathways,
    pathways
  ) %>%
  slice_head(n = 10)

top10_combined <- bind_rows(
  top10_gwas_genes %>% mutate(table_group = "Top10_GWAS_genes", .before = 1),
  top10_magma_genes %>% mutate(table_group = "Top10_MAGMA_genes", .before = 1),
  top10_all3_genes %>% mutate(table_group = "Top10_All3_genes", .before = 1)
)

write_me(top10_gwas_genes, "SuppTable_DryTons_top10_GWAS_genes_from_scoring.csv")
write_me(top10_magma_genes, "SuppTable_DryTons_top10_MAGMA_genes_from_scoring.csv")
write_me(top10_all3_genes, "SuppTable_DryTons_top10_All3_candidate_genes_from_scoring.csv")
write_me(top10_combined, "SuppTable_DryTons_top10_gene_panels_from_scoring_combined.csv")

fig7_candidate_summary <- candidate_with_class %>%
  count(evidence_class, name = "n_genes") %>%
  arrange(desc(n_genes), evidence_class)

write_me(fig7_candidate_summary, "SuppTable_Fig7_candidate_evidence_summary_Dry_tons_per_acre.csv")

# ------------------------------------------------------------------------------
# Lignin: top 50 candidate genes
# ------------------------------------------------------------------------------
lignin_candidate_file <- file.path(
  "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Lignin_DM_perc/candidate_gene_scoring_B10000_GPD_strict_tau",
  "candidate_genes_top50_Lignin_DM_perc_B10000_GPD_strict_tau.csv"
)

if (file.exists(lignin_candidate_file)) {
  lignin_top50 <- read_csv(lignin_candidate_file, show_col_types = FALSE)
  write_me(lignin_top50, "SuppTable_Lignin_top50_candidate_genes.csv")
}

# ------------------------------------------------------------------------------
# README
# ------------------------------------------------------------------------------
readme_lines <- c(
  "Supplementary tables linked to final manuscript figures.",
  "",
  "Figures covered:",
  "- Fig1: simulated archetype rank curves",
  "- Fig3: pathway-size sensitivity / archetype verification",
  "- Fig4: null calibration and missing-correlation calibration",
  "- Fig5: component breakage and adaptive / leave-one-out diagnostics",
  "- Fig6: Dry tons per acre GWAS, MAGMA, and CATFISH pathway integration",
  "- Fig7: Dry tons per acre CATFISH component comparison and candidate prioritization",
  "- Extra trait table: Lignin top 50 candidate genes"
)
writeLines(readme_lines, con = file.path(out_dir, "README.txt"))

cat("Supplementary tables written to:\n", out_dir, "\n", sep = "")
