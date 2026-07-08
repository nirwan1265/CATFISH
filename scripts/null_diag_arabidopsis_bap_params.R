#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x) || !nzchar(x)) y else x

suppressPackageStartupMessages({
  library(devtools)
  library(CATFISH)
  library(dplyr)
  library(data.table)
})

this_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
this_file <- if (length(this_arg)) sub("^--file=", "", this_arg[[1]]) else file.path(getwd(), "scripts/null_diag_arabidopsis_bap_params.R")
repo_dir <- normalizePath(file.path(dirname(this_file), ".."), mustWork = TRUE)
devtools::load_all(repo_dir, quiet = TRUE)

get_env_num <- function(name, default) {
  x <- Sys.getenv(name, "")
  if (!nzchar(x)) return(default)
  out <- suppressWarnings(as.numeric(x))
  if (!is.finite(out) || is.na(out)) default else out
}

ARAB_MAGMA_DIR <- Sys.getenv(
  "ARAB_MAGMA_DIR",
  "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/MAGMA/AT_cold_by_chr"
)
ARAB_COR_FILE <- Sys.getenv(
  "ARAB_COR_FILE",
  "/Users/nirwantandukar/Documents/Research/results/MAGMA/MAGCAT/magma_multi_snp_wise_genes_by_chr_N_maize/magma_gene_cor_pairs_MLM_arabidopsis.txt"
)
ARAB_GENE_LEN <- Sys.getenv(
  "ARAB_GENE_LEN",
  file.path(repo_dir, "inst/extdata/Arabidopsis_gene_lengths.tsv")
)
ARAB_PATHWAY_FILE <- Sys.getenv(
  "ARAB_PATHWAY_FILE",
  file.path(repo_dir, "inst/extdata/pathway/aracyc_pathways.20230103")
)
OUT_DIR <- Sys.getenv(
  "ARAB_NULL_OUT",
  file.path(repo_dir, "null_diagnostics", "arabidopsis_bap_params")
)

N_NULL <- as.integer(get_env_num("N_NULL", 20))
B_PERM <- as.integer(get_env_num("B_PERM", 10000))
SEED <- as.integer(get_env_num("SEED", 42))
TAIL_GPD_K <- as.integer(get_env_num("TAIL_GPD_K", 250))
TAIL_MIN_B <- as.integer(get_env_num("TAIL_MIN_B", 1000))
N_THREADS <- as.integer(get_env_num("CATFISH_THREADS", 1))

TAU_GRID <- c(0.10, 0.05, 0.01)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

compute_lambda <- function(p_values) {
  p <- p_values[is.finite(p_values) & !is.na(p_values) & p_values > 0 & p_values < 1]
  if (length(p) < 10L) return(NA_real_)
  median(qchisq(1 - p, df = 1)) / qchisq(0.5, df = 1)
}

compute_type1_error <- function(p_values, alpha = 0.05) {
  p <- p_values[is.finite(p_values) & !is.na(p_values)]
  if (!length(p)) return(NA_real_)
  mean(p < alpha)
}

generate_null_genes <- function(gene_template, min_p = 1e-10, max_p = 1 - 1e-10) {
  out <- gene_template
  out$P <- runif(nrow(out))
  out$P <- pmax(pmin(out$P, max_p), min_p)
  out$ZSTAT <- qnorm(1 - out$P)
  out$ZSTAT[!is.finite(out$ZSTAT)] <- 0
  out
}

cat("[arab-null] Loading MAGMA genes from:", ARAB_MAGMA_DIR, "\n")
files <- list.files(
  path = ARAB_MAGMA_DIR,
  pattern = "^AT_cold_chr.*\\.multi_snp_wise\\.genes\\.out$",
  full.names = TRUE
)
if (!length(files)) stop("No Arabidopsis MAGMA gene files found in ", ARAB_MAGMA_DIR, call. = FALSE)

gene_list <- lapply(files, function(f) {
  read.table(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
})
genes_all_raw <- do.call(rbind, gene_list)
colnames(genes_all_raw)[9] <- "P"
o <- order(genes_all_raw$GENE, genes_all_raw$P)
genes_all <- genes_all_raw[o, ]
genes_all <- genes_all[!duplicated(genes_all$GENE), ]

gene_len <- read.delim(ARAB_GENE_LEN, stringsAsFactors = FALSE)
gene_len$gene_id <- sub("^gene:", "", gene_len$gene_id)

adj_out <- catfish_adjust_gene_p(
  gene_results = genes_all,
  gene_lengths = gene_len,
  gene_col = "GENE",
  nsnp_col = "NSNPS",
  p_col = "P",
  z_col = "ZSTAT",
  len_gene_col = "gene_id",
  len_col = "length"
)

gene_template <- data.frame(
  GENE = adj_out[[1]],
  ZSTAT = adj_out[[2]],
  P = adj_out[[3]],
  stringsAsFactors = FALSE
)

pathways <- read.delim(ARAB_PATHWAY_FILE, stringsAsFactors = FALSE) %>%
  transmute(
    pathway_id = Pathway.id,
    pathway_name = Pathway.name,
    gene_id = sub("^gene:", "", Gene.id)
  ) %>%
  filter(!is.na(gene_id), gene_id != "") %>%
  distinct()

cat("[arab-null] Settings:\n")
cat("  N_NULL =", N_NULL, "\n")
cat("  B_PERM =", B_PERM, "\n")
cat("  TAU_GRID =", paste(TAU_GRID, collapse = ","), "\n")
cat("  mvn_calibrate_components = FALSE\n")
cat("  tail_mode = hybrid_gpd\n")
cat("  n_threads =", N_THREADS, "\n")

keep_cols <- c(
  "acat_p",
  "fisher_p",
  "tfisher_p_analytic",
  "minp_p_analytic",
  "stouffer_p_analytic",
  "omni_p_analytic",
  "omni_p_mvn",
  "omni_p_final"
)

results <- vector("list", N_NULL)
for (i in seq_len(N_NULL)) {
  cat(sprintf("[arab-null] iter %d/%d\n", i, N_NULL))
  null_genes <- generate_null_genes(gene_template)
  res <- catfish_omni2_pathways(
    gene_results = null_genes,
    species = NULL,
    pathways = pathways,
    pmn_gene_col = "Gene-name",
    gene_col = "GENE",
    p_raw_col = "P",
    z_col = "ZSTAT",
    tau_grid = TAU_GRID,
    min_p = 1e-15,
    do_fix = TRUE,
    stouffer_min_abs_w = 1e-8,
    stouffer_alternative = "greater",
    include_magma_in_omni = FALSE,
    include_magma_in_perm = FALSE,
    omnibus = "ACAT",
    B_perm = B_PERM,
    perm_mode = "mvn_global",
    magma_cor_file = ARAB_COR_FILE,
    make_PD = TRUE,
    mvn_marginal = "uniform",
    mvn_calibrate_components = FALSE,
    tail_mode = "hybrid_gpd",
    tail_gpd_k = TAIL_GPD_K,
    tail_min_B = TAIL_MIN_B,
    output = FALSE,
    seed = SEED + i,
    n_threads = N_THREADS
  )
  results[[i]] <- res[, intersect(c("pathway_id", keep_cols), names(res))]
  results[[i]]$null_iter <- i
}

null_df <- bind_rows(results)
fwrite(null_df, file.path(OUT_DIR, "null_pvalues_arabidopsis_bap_params.csv"))

summary_rows <- lapply(intersect(keep_cols, names(null_df)), function(col) {
  p <- null_df[[col]]
  data.frame(
    col = col,
    n_pvals = sum(!is.na(p)),
    lambda = compute_lambda(p),
    type1_alpha005 = compute_type1_error(p, 0.05),
    type1_alpha001 = compute_type1_error(p, 0.01),
    mean_p = mean(p, na.rm = TRUE),
    median_p = median(p, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
})
summary_df <- bind_rows(summary_rows)
fwrite(summary_df, file.path(OUT_DIR, "summary_arabidopsis_bap_params.csv"))

cat("[arab-null] Wrote:\n")
cat("  ", file.path(OUT_DIR, "null_pvalues_arabidopsis_bap_params.csv"), "\n")
cat("  ", file.path(OUT_DIR, "summary_arabidopsis_bap_params.csv"), "\n")
print(summary_df)
