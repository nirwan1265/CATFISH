#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(readr)
})

## -----------------------------------------------------------------------------
## Parameters
## -----------------------------------------------------------------------------

GWAS_MODE           <- "top_pct"
GWAS_P_THRESHOLD    <- 1e-5
GWAS_TOP_PCT        <- 1

MAGMA_MODE          <- "top_pct"
MAGMA_FDR_THRESHOLD <- 0.30
MAGMA_TOP_PCT       <- 5

WINDOW_SIZE <- 25000

BASE_DIR <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre"
GWAS_DIR <- "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/BAP/Dry_tons_per_acre"
GENE_LOC_FILE <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/sorghum.genes.loc"

CATFISH_SUBDIR_OVERRIDE <- Sys.getenv("CATFISH_SUBDIR", "")
CATFISH_CSV_OVERRIDE <- Sys.getenv("CATFISH_CSV", "")
MAGMA_GENE_FILE_OVERRIDE <- Sys.getenv("MAGMA_GENE_FILE", "")
OUT_LABEL <- Sys.getenv("OUT_LABEL", "B1000000_GPD_strict_tau")

OUT_DIR <- file.path(BASE_DIR, paste0("candidate_gene_scoring_", OUT_LABEL))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------

pick_completed_catfish_csv <- function(base_dir) {
  if (nzchar(CATFISH_CSV_OVERRIDE)) {
    if (!file.exists(CATFISH_CSV_OVERRIDE)) {
      stop("CATFISH_CSV override not found: ", CATFISH_CSV_OVERRIDE, call. = FALSE)
    }
    return(CATFISH_CSV_OVERRIDE)
  }
  if (nzchar(CATFISH_SUBDIR_OVERRIDE)) {
    hits <- list.files(
      file.path(base_dir, CATFISH_SUBDIR_OVERRIDE),
      pattern = "^Dry_tons_per_acre_CATFISH_.*\\.csv$",
      full.names = TRUE
    )
    if (!length(hits)) {
      stop("No CATFISH CSV found in subdir override: ", CATFISH_SUBDIR_OVERRIDE, call. = FALSE)
    }
    return(hits[[1]])
  }
  csvs <- c(
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn_GPD_strict_tau", "Dry_tons_per_acre_CATFISH_ACAT_mvn_B1000000_GPD_strict_tau.csv"),
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn_GPD", "Dry_tons_per_acre_CATFISH_ACAT_mvn_B1000000_GPD.csv"),
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn", "Dry_tons_per_acre_CATFISH_ACAT_mvn_B1000000.csv"),
    file.path(base_dir, "CATFISH_permutation_B10000_mvn", "Dry_tons_per_acre_CATFISH_ACAT_mvn_B10000.csv")
  )
  hit <- csvs[file.exists(csvs)]
  if (!length(hit)) {
    stop("No completed CATFISH result CSV found under: ", base_dir, call. = FALSE)
  }
  hit[[1]]
}

pick_completed_gene_table <- function(base_dir) {
  if (nzchar(MAGMA_GENE_FILE_OVERRIDE)) {
    if (!file.exists(MAGMA_GENE_FILE_OVERRIDE)) {
      stop("MAGMA_GENE_FILE override not found: ", MAGMA_GENE_FILE_OVERRIDE, call. = FALSE)
    }
    return(MAGMA_GENE_FILE_OVERRIDE)
  }
  if (nzchar(CATFISH_SUBDIR_OVERRIDE)) {
    hit <- file.path(base_dir, CATFISH_SUBDIR_OVERRIDE, "Dry_tons_per_acre_combined_genes.tsv")
    if (!file.exists(hit)) {
      stop("Combined gene table not found in subdir override: ", CATFISH_SUBDIR_OVERRIDE, call. = FALSE)
    }
    return(hit)
  }
  tsvs <- c(
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn_GPD_strict_tau", "Dry_tons_per_acre_combined_genes.tsv"),
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn_GPD", "Dry_tons_per_acre_combined_genes.tsv"),
    file.path(base_dir, "CATFISH_permutation_B1000000_mvn", "Dry_tons_per_acre_combined_genes.tsv"),
    file.path(base_dir, "CATFISH_permutation_B10000_mvn", "Dry_tons_per_acre_combined_genes.tsv")
  )
  hit <- tsvs[file.exists(tsvs)]
  if (!length(hit)) {
    stop("No combined MAGMA gene table found under: ", base_dir, call. = FALSE)
  }
  hit[[1]]
}

read_gene_loc_clean <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("GENE", "CHR", "START", "STOP")
  miss <- setdiff(req, names(x))
  if (length(miss)) {
    stop("Gene location file missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }
  x <- x[, req, drop = FALSE]
  x$GENE <- as.character(x$GENE)
  x$CHR <- as.integer(x$CHR)
  x$START <- as.integer(x$START)
  x$STOP <- as.integer(x$STOP)
  x <- x[!duplicated(x$GENE), , drop = FALSE]
  x$START_EXT <- pmax(0L, x$START - WINDOW_SIZE)
  x$STOP_EXT  <- x$STOP + WINDOW_SIZE
  x
}

read_all_gwas <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (!length(files)) stop("No GWAS .assoc.txt files found in: ", dir_path, call. = FALSE)

  chr_num <- as.integer(sub(".*Chr([0-9]+)\\.maf01\\.assoc\\.txt$", "\\1", basename(files)))
  files <- files[order(chr_num)]

  pieces <- lapply(files, function(f) {
    x <- fread(f, header = TRUE)
    names(x)[names(x) == "chr"] <- "CHR"
    names(x)[names(x) == "rs"] <- "SNP_ID"
    names(x)[names(x) == "ps"] <- "POS"
    names(x)[names(x) == "p_wald"] <- "P"
    names(x)[names(x) == "n_miss"] <- "NMISS"
    x
  })

  out <- rbindlist(pieces, use.names = TRUE, fill = TRUE)
  req <- c("CHR", "SNP_ID", "POS", "P")
  miss <- setdiff(req, names(out))
  if (length(miss)) {
    stop("GWAS files missing required columns: ", paste(miss, collapse = ", "), call. = FALSE)
  }

  out[, CHR := as.integer(CHR)]
  out[, POS := as.integer(POS)]
  out
}

## -----------------------------------------------------------------------------
## Inputs
## -----------------------------------------------------------------------------

CATFISH_CSV <- pick_completed_catfish_csv(BASE_DIR)
MAGMA_GENE_FILE <- pick_completed_gene_table(BASE_DIR)

cat("Using CATFISH CSV: ", CATFISH_CSV, "\n", sep = "")
cat("Using MAGMA gene table: ", MAGMA_GENE_FILE, "\n", sep = "")
cat("GWAS dir: ", GWAS_DIR, "\n", sep = "")
cat("Output dir: ", OUT_DIR, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 1. Gene annotation
## -----------------------------------------------------------------------------

gene_loc <- read_gene_loc_clean(GENE_LOC_FILE)
cat("Genes in annotation after dedup: ", nrow(gene_loc), "\n", sep = "")

## -----------------------------------------------------------------------------
## 2. GWAS layer
## -----------------------------------------------------------------------------

gwas_dt <- read_all_gwas(GWAS_DIR)
cat("Total GWAS SNPs: ", nrow(gwas_dt), "\n", sep = "")

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

gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_n_snps = n(),
    .groups = "drop"
  ) %>%
  arrange(gwas_min_p, GENE) %>%
  mutate(gwas_rank = dplyr::row_number())

if (GWAS_MODE == "top_pct") {
  gwas_n_top <- ceiling(nrow(gwas_gene) * GWAS_TOP_PCT / 100)
  gwas_eff_threshold <- sort(gwas_gene$gwas_min_p)[gwas_n_top]
} else {
  gwas_eff_threshold <- GWAS_P_THRESHOLD
}

cat("Genes with GWAS evidence: ", nrow(gwas_gene), "\n", sep = "")
cat("GWAS effective threshold: ", format(gwas_eff_threshold, digits = 4), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 3. MAGMA layer
## -----------------------------------------------------------------------------

magma_gene <- read.delim(MAGMA_GENE_FILE, stringsAsFactors = FALSE, check.names = FALSE)
if (!"P" %in% names(magma_gene) && "P_MULTI" %in% names(magma_gene)) magma_gene$P <- magma_gene$P_MULTI
if (!all(c("GENE", "P", "ZSTAT") %in% names(magma_gene))) {
  stop("MAGMA gene table must contain GENE, P, ZSTAT.", call. = FALSE)
}

magma_gene <- magma_gene %>%
  rename(magma_p = P, magma_z = ZSTAT) %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_bonf = p.adjust(magma_p, method = "bonferroni")
  ) %>%
  arrange(magma_p, GENE) %>%
  mutate(magma_rank = dplyr::row_number())

if (MAGMA_MODE == "top_pct") {
  magma_n_top <- ceiling(nrow(magma_gene) * MAGMA_TOP_PCT / 100)
  magma_eff_threshold <- sort(magma_gene$magma_p)[magma_n_top]
} else {
  magma_eff_threshold <- NULL
}

cat("Genes with MAGMA results: ", nrow(magma_gene), "\n", sep = "")
if (MAGMA_MODE == "top_pct") {
  cat("MAGMA effective threshold: ", format(magma_eff_threshold, digits = 4), "\n\n", sep = "")
}

## -----------------------------------------------------------------------------
## 4. Pathway layer
## -----------------------------------------------------------------------------

omni_results <- read.csv(CATFISH_CSV, stringsAsFactors = FALSE)

omni_col <- if ("omni_p_final" %in% names(omni_results)) {
  "omni_p_final"
} else if ("omni_p_mvn" %in% names(omni_results)) {
  "omni_p_mvn"
} else {
  "omni_p_analytic"
}

omni_results <- omni_results %>%
  arrange(.data[[omni_col]], pathway_name, pathway_id) %>%
  mutate(pathway_rank = dplyr::row_number())

pathway_genes_list <- strsplit(omni_results$genes_used, ";")
pathway_genes <- data.frame(
  pathway_id = rep(omni_results$pathway_id, lengths(pathway_genes_list)),
  pathway_name = rep(omni_results$pathway_name, lengths(pathway_genes_list)),
  pathway_p = rep(omni_results[[omni_col]], lengths(pathway_genes_list)),
  pathway_rank = rep(omni_results$pathway_rank, lengths(pathway_genes_list)),
  GENE = trimws(unlist(pathway_genes_list)),
  stringsAsFactors = FALSE
)

pathway_genes <- pathway_genes %>%
  filter(!is.na(GENE), nzchar(GENE))

pathway_best <- pathway_genes %>%
  arrange(pathway_rank, pathway_name, pathway_id) %>%
  group_by(GENE) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    GENE,
    best_pathway_id = pathway_id,
    best_pathway_name = pathway_name,
    best_pathway_p = pathway_p,
    pathway_rank = pathway_rank
  )

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(
    n_top_pathways = n_distinct(pathway_id),
    mean_pathway_mlog10p = mean(-log10(pathway_p), na.rm = TRUE),
    pathways = paste(pathway_id, collapse = "; "),
    pathway_names = paste(pathway_name, collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(pathway_best, by = "GENE") %>%
  arrange(pathway_rank, GENE)

n_pathways_total <- nrow(omni_results)
cat("Total pathways used for scoring: ", n_pathways_total, "\n", sep = "")
cat("Genes with pathway membership: ", nrow(pathway_gene_support), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 5. Integrate three layers and score genes
## -----------------------------------------------------------------------------

n_genes_total <- nrow(magma_gene)

gene_evidence <- magma_gene %>%
  left_join(gwas_gene, by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas = !is.na(gwas_min_p) & gwas_min_p <= gwas_eff_threshold,
    hit_magma = if (MAGMA_MODE == "top_pct") {
      !is.na(magma_p) & magma_p <= magma_eff_threshold
    } else {
      !is.na(magma_fdr) & magma_fdr < MAGMA_FDR_THRESHOLD
    },
    hit_pathway = !is.na(n_top_pathways) & n_top_pathways >= 1,
    support_layers = hit_gwas + hit_magma + hit_pathway,
    gwas_rank_score = ifelse(!is.na(gwas_rank), -log10(gwas_rank / n_genes_total), 0),
    magma_rank_score = ifelse(!is.na(magma_rank), -log10(magma_rank / n_genes_total), 0),
    pathway_rank_score = ifelse(!is.na(pathway_rank), -log10(pathway_rank / n_pathways_total), 0),
    multi_pathway_bonus = ifelse(!is.na(n_top_pathways), log10(n_top_pathways + 1), 0),
    score = gwas_rank_score +
      magma_rank_score +
      0.5 * pathway_rank_score +
      0.3 * multi_pathway_bonus
  ) %>%
  arrange(desc(score), gwas_rank, magma_rank, pathway_rank, GENE)

top50 <- gene_evidence %>%
  select(
    GENE, score, support_layers,
    gwas_rank_score, magma_rank_score, pathway_rank_score, multi_pathway_bonus,
    magma_p, magma_fdr, magma_rank, magma_z,
    gwas_min_p, gwas_n_snps, gwas_rank,
    n_top_pathways, best_pathway_p, pathway_rank, best_pathway_id, best_pathway_name, pathways, pathway_names,
    hit_gwas, hit_magma, hit_pathway
  ) %>%
  slice_head(n = 50)

top200 <- gene_evidence %>%
  select(
    GENE, score, support_layers,
    gwas_rank_score, magma_rank_score, pathway_rank_score, multi_pathway_bonus,
    magma_p, magma_fdr, magma_rank, magma_z,
    gwas_min_p, gwas_n_snps, gwas_rank,
    n_top_pathways, best_pathway_p, pathway_rank, best_pathway_id, best_pathway_name, pathways, pathway_names,
    hit_gwas, hit_magma, hit_pathway
  ) %>%
  slice_head(n = 200)

all_file <- file.path(OUT_DIR, paste0("candidate_genes_all_Dry_tons_per_acre_", OUT_LABEL, ".csv"))
top50_file <- file.path(OUT_DIR, paste0("candidate_genes_top50_Dry_tons_per_acre_", OUT_LABEL, ".csv"))
top200_file <- file.path(OUT_DIR, paste0("candidate_genes_top200_Dry_tons_per_acre_", OUT_LABEL, ".csv"))

write.csv(gene_evidence, all_file, row.names = FALSE)
write.csv(top50, top50_file, row.names = FALSE)
write.csv(top200, top200_file, row.names = FALSE)

cat("Saved files:\n")
cat("  ", all_file, "\n", sep = "")
cat("  ", top50_file, "\n", sep = "")
cat("  ", top200_file, "\n\n", sep = "")

cat("Top 10 genes by score:\n")
print(top50 %>% slice_head(n = 10))
