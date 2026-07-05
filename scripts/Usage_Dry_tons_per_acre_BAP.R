#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
})
load_all(".")

default_repo <- Sys.getenv("CATFISH_REPO", "/Users/nirwantandukar/Documents/Github/MAGCAT")

## -----------------------------------------------------------------------------
## Configuration
## -----------------------------------------------------------------------------

get_env_int <- function(name, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, "")))
  if (!is.finite(x) || is.na(x)) default else x
}

get_env_flag <- function(name, default = FALSE) {
  raw <- Sys.getenv(name, "")
  if (!nzchar(raw)) return(default)
  tolower(trimws(raw)) %in% c("1", "true", "yes", "y")
}

parse_num_grid <- function(x, default) {
  if (!nzchar(x)) return(default)
  parts <- unlist(strsplit(gsub("[,;]", " ", x), "\\s+"))
  vals <- suppressWarnings(as.numeric(parts))
  vals <- vals[is.finite(vals) & !is.na(vals)]
  if (!length(vals)) default else vals
}

MAGMA_DIR <- Sys.getenv(
  "MAGMA_DIR",
  "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre"
)
PATHWAY_FILE <- Sys.getenv(
  "PATHWAY_FILE",
  file.path(default_repo, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI")
)

OUT_PREFIX <- "Dry_tons_per_acre"
SPECIES    <- NULL
GENE_REGEX <- "^SORBI"
TAU_OPTION <- Sys.getenv("TAU_OPTION", "paper")

B_PERM      <- get_env_int("B_PERM", 1000000L)
PERM_MODE   <- Sys.getenv("PERM_MODE", "mvn")
OMNIBUS     <- Sys.getenv("OMNIBUS", "ACAT")
SEED        <- get_env_int("SEED", 123L)
N_THREADS   <- get_env_int("N_THREADS", min(12L, max(1L, parallel::detectCores() - 1L)))
MIN_GENES   <- get_env_int("MIN_GENES", 2L)

tau_grid_default <- switch(
  TAU_OPTION,
  strict  = c(1e-5, 1e-6, 1e-7),
  paper   = c(0.1, 0.05, 0.01),
  default = c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001),
  c(0.1, 0.05, 0.02, 0.01, 0.005, 0.001)
)
TAU_GRID <- parse_num_grid(Sys.getenv("TAU_GRID", ""), tau_grid_default)
TAU_LABEL <- Sys.getenv(
  "TAU_LABEL",
  if (identical(TAU_OPTION, "strict")) "strict_tau" else if (identical(TAU_OPTION, "paper")) "paper_tau" else "default_tau"
)
OUT_SUFFIX <- if (nzchar(TAU_LABEL)) paste0("_", TAU_LABEL) else ""

MVN_MARGINAL             <- Sys.getenv("MVN_MARGINAL", "uniform")
MVN_CALIBRATE_COMPONENTS <- get_env_flag("MVN_CALIBRATE_COMPONENTS", TRUE)
MAKE_PD                  <- get_env_flag("MAKE_PD", TRUE)
TAIL_MODE                <- Sys.getenv("TAIL_MODE", "hybrid_gpd")
TAIL_SWITCH_EXCEED       <- get_env_int("TAIL_SWITCH_EXCEED", 10L)
TAIL_GPD_K               <- get_env_int("TAIL_GPD_K", 250L)
TAIL_MIN_B               <- get_env_int("TAIL_MIN_B", 10000L)
TAIL_MIN_TAIL            <- get_env_int("TAIL_MIN_TAIL", 50L)

tail_label <- if (identical(TAIL_MODE, "hybrid_gpd")) "GPD" else "empirical"
OUT_DIR <- file.path(MAGMA_DIR, paste0("CATFISH_permutation_B", B_PERM, "_", PERM_MODE, "_", tail_label, OUT_SUFFIX))

USE_ADJUSTED_GENE_P <- FALSE
GENE_LENGTH_FILE    <- NA_character_

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------

read_magma_gene_out <- function(path) {
  x <- utils::read.table(
    path,
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  if (!"GENE" %in% names(x)) {
    stop("Missing GENE column in: ", path, call. = FALSE)
  }
  if (!"ZSTAT" %in% names(x)) {
    stop("Missing ZSTAT column in: ", path, call. = FALSE)
  }

  p_col <- if ("P_MULTI" %in% names(x)) "P_MULTI" else if ("P" %in% names(x)) "P" else NA_character_
  if (is.na(p_col)) {
    stop("Could not find a MAGMA p-value column in: ", path, call. = FALSE)
  }

  x$P <- x[[p_col]]
  x
}

combine_gene_results <- function(files) {
  gene_list <- lapply(files, read_magma_gene_out)
  genes_all <- do.call(rbind, gene_list)
  genes_all <- genes_all[order(genes_all$GENE, genes_all$P), , drop = FALSE]
  genes_all <- genes_all[!duplicated(genes_all$GENE), , drop = FALSE]
  rownames(genes_all) <- NULL
  genes_all
}

write_combined_cor_pairs <- function(raw_files, out_file, gene_regex) {
  if (file.exists(out_file)) file.remove(out_file)

  first <- TRUE
  n_written <- 0L

  for (f in raw_files) {
    tmp <- tempfile(fileext = ".txt")

    magma_genesraw_to_cor_pairs_banded(
      genes_raw_file = f,
      out_pairs_file = tmp,
      gene_regex     = gene_regex,
      keep_abs_r_ge  = 0,
      overwrite      = TRUE,
      verbose        = FALSE
    )

    x <- readLines(tmp, warn = FALSE)
    x <- x[nzchar(x)]
    if (!length(x)) next
    if (!first) x <- x[-1]
    if (!length(x)) next

    cat(paste(x, collapse = "\n"), "\n", file = out_file, append = !first, sep = "")
    n_written <- n_written + length(x)
    first <- FALSE
  }

  if (!file.exists(out_file)) {
    stop("Failed to create correlation-pair file: ", out_file, call. = FALSE)
  }

  invisible(n_written)
}

pick_final_p_col <- function(x) {
  picks <- c("omni_p_final", "omni_p_mvn", "omni_p_global", "omni_p_analytic")
  hit <- picks[picks %in% names(x)]
  if (length(hit)) hit[[1]] else NA_character_
}

## -----------------------------------------------------------------------------
## Inputs
## -----------------------------------------------------------------------------

gene_out_files <- list.files(
  MAGMA_DIR,
  pattern = paste0("^", OUT_PREFIX, "_chr[0-9]+\\.multi_snp_wise\\.genes\\.out$"),
  full.names = TRUE
)

gene_raw_files <- list.files(
  MAGMA_DIR,
  pattern = paste0("^", OUT_PREFIX, "_chr[0-9]+\\.multi_snp_wise\\.genes\\.raw$"),
  full.names = TRUE
)

if (!length(gene_out_files)) {
  stop("No per-chromosome .genes.out files found in: ", MAGMA_DIR, call. = FALSE)
}
if (!length(gene_raw_files)) {
  stop("No per-chromosome .genes.raw files found in: ", MAGMA_DIR, call. = FALSE)
}

chr_num_out <- as.integer(sub(".*_chr([0-9]+)\\..*$", "\\1", basename(gene_out_files)))
gene_out_files <- gene_out_files[order(chr_num_out)]

chr_num_raw <- as.integer(sub(".*_chr([0-9]+)\\..*$", "\\1", basename(gene_raw_files)))
gene_raw_files <- gene_raw_files[order(chr_num_raw)]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("MAGMA dir: ", MAGMA_DIR, "\n", sep = "")
cat("Output dir: ", OUT_DIR, "\n", sep = "")
cat("Gene .out files: ", length(gene_out_files), "\n", sep = "")
cat("Gene .raw files: ", length(gene_raw_files), "\n", sep = "")
cat("Permutation mode: ", PERM_MODE, "\n", sep = "")
cat("B_perm: ", B_PERM, "\n", sep = "")
cat("Requested cores for surrounding work: ", N_THREADS, "\n\n", sep = "")
cat("Tail mode: ", TAIL_MODE, "\n", sep = "")
cat("MVN calibrate components: ", MVN_CALIBRATE_COMPONENTS, "\n\n", sep = "")
cat("Tau option: ", TAU_OPTION, "\n", sep = "")
cat("Tau grid: ", paste(format(TAU_GRID, scientific = TRUE), collapse = ", "), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## Step 1: Combine gene-level MAGMA results
## -----------------------------------------------------------------------------

genes_all <- combine_gene_results(gene_out_files)

combined_gene_file <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_combined_genes.tsv"))
utils::write.table(
  genes_all,
  combined_gene_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Combined genes written to: ", combined_gene_file, "\n", sep = "")
cat("Unique genes: ", nrow(genes_all), "\n\n", sep = "")

genes_for_catfish <- genes_all

if (USE_ADJUSTED_GENE_P) {
  if (!file.exists(GENE_LENGTH_FILE)) {
    stop("GENE_LENGTH_FILE not found: ", GENE_LENGTH_FILE, call. = FALSE)
  }

  gene_lengths <- utils::read.table(
    GENE_LENGTH_FILE,
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  adj_out <- CATFISH::catfish_adjust_gene_p(
    gene_results = genes_all,
    gene_lengths = gene_lengths,
    gene_col     = "GENE",
    nsnp_col     = "NSNPS",
    p_col        = "P",
    z_col        = "ZSTAT",
    len_gene_col = "gene_id",
    len_col      = "length"
  )

  genes_for_catfish <- genes_all
  genes_for_catfish$P_adj <- adj_out$p_adj[match(genes_all$GENE, adj_out$gene_id)]
}

## -----------------------------------------------------------------------------
## Step 2: Build gene-gene correlation pairs from MAGMA .genes.raw
## -----------------------------------------------------------------------------

cor_pairs_file <- file.path(OUT_DIR, paste0(OUT_PREFIX, "_gene_cor_pairs.txt"))
n_pairs_lines <- write_combined_cor_pairs(
  raw_files  = gene_raw_files,
  out_file   = cor_pairs_file,
  gene_regex = GENE_REGEX
)

cat("Correlation-pair file written to: ", cor_pairs_file, "\n", sep = "")
cat("Correlation lines written: ", n_pairs_lines, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## Step 3: Load SORBI-format sorghum pathways
## -----------------------------------------------------------------------------

if (!file.exists(PATHWAY_FILE)) {
  stop("PATHWAY_FILE not found: ", PATHWAY_FILE, call. = FALSE)
}

pathway_raw <- utils::read.delim(PATHWAY_FILE, header = TRUE, stringsAsFactors = FALSE)
pathways <- unique(data.frame(
  pathway_id   = pathway_raw[["Pathway.id"]],
  pathway_name = pathway_raw[["Pathway.name"]],
  gene_id      = pathway_raw[["Gene.name"]],
  stringsAsFactors = FALSE
))

cat("Loaded pathways: ", length(unique(pathways$pathway_id)), "\n", sep = "")
cat("Pathway-gene rows: ", nrow(pathways), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## Step 4: Run CATFISH with MVN permutation calibration
## -----------------------------------------------------------------------------

omni_results <- catfish_omni2_pathways(
  gene_results              = genes_for_catfish,
  species                   = SPECIES,
  pathways                  = pathways,
  gene_col                  = "GENE",
  p_raw_col                 = "P",
  p_adj_col                 = if (USE_ADJUSTED_GENE_P) "P_adj" else NULL,
  z_col                     = "ZSTAT",
  tau_grid                  = TAU_GRID,
  min_p                     = 1e-15,
  do_fix                    = TRUE,
  stouffer_alternative      = "greater",
  omnibus                   = OMNIBUS,
  perm_mode                 = PERM_MODE,
  B_perm                    = B_PERM,
  magma_cor_file            = cor_pairs_file,
  make_PD                   = MAKE_PD,
  mvn_marginal              = MVN_MARGINAL,
  mvn_calibrate_components  = MVN_CALIBRATE_COMPONENTS,
  tail_mode                 = TAIL_MODE,
  tail_switch_exceed        = TAIL_SWITCH_EXCEED,
  tail_gpd_k                = TAIL_GPD_K,
  tail_min_B                = TAIL_MIN_B,
  tail_min_tail             = TAIL_MIN_TAIL,
  n_threads                 = N_THREADS,
  min_genes                 = MIN_GENES,
  seed                      = SEED,
  output                    = FALSE
)

final_p_col <- pick_final_p_col(omni_results)
if (is.na(final_p_col)) {
  stop("Could not find a final omnibus p-value column in CATFISH output.", call. = FALSE)
}

omni_results$final_p_col_used <- final_p_col
omni_results$FDR_BH <- stats::p.adjust(omni_results[[final_p_col]], method = "BH")
omni_results <- omni_results[order(omni_results[[final_p_col]], omni_results$FDR_BH), , drop = FALSE]

results_csv <- file.path(
  OUT_DIR,
  paste0(OUT_PREFIX, "_CATFISH_", OMNIBUS, "_", PERM_MODE, "_B", B_PERM, "_", tail_label, OUT_SUFFIX, ".csv")
)
results_rds <- file.path(
  OUT_DIR,
  paste0(OUT_PREFIX, "_CATFISH_", OMNIBUS, "_", PERM_MODE, "_B", B_PERM, "_", tail_label, OUT_SUFFIX, ".rds")
)

utils::write.csv(omni_results, results_csv, row.names = FALSE)
saveRDS(omni_results, results_rds)

cat("CATFISH results written to:\n")
cat("  ", results_csv, "\n", sep = "")
cat("  ", results_rds, "\n\n", sep = "")

cat("Top pathways:\n")
print(utils::head(omni_results[, intersect(
  c("pathway_id", "pathway_name", "omni_p_final", "omni_p_mvn", "omni_p_analytic", "FDR_BH"),
  names(omni_results)
)], 10))
