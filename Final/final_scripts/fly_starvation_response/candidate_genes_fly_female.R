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

TRAIT    <- "Starvation_stress_female"
BASE_DIR <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/DGRP/Fly_starvation_female"
GWAS_FILE <- "/Users/nirwantandukar/Documents/Research/data/DGRP/Starvation_stress/raw_gwas/raw_GWAS_Starvation_stress_female_DGRP.csv"
GENE_LOC_FILE <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/fly.genes.loc"

OUT_DIR <- file.path(BASE_DIR, "candidate_gene_scoring_B100000_GPD")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Chromosome arm → numeric code (matches MAGMA gene loc encoding)
CHR_MAP <- c("2L" = 1L, "2R" = 2L, "3L" = 3L, "3R" = 4L, "4" = 5L, "X" = 6L, "Y" = 7L)

## -----------------------------------------------------------------------------
## Helpers
## -----------------------------------------------------------------------------

pick_completed_catfish_csv <- function(base_dir) {
  csvs <- c(
    file.path(base_dir, "CATFISH_permutation_B100000_mvn_GPD",
              paste0(TRAIT, "_CATFISH_ACAT_mvn_B100000_GPD.csv")),
    file.path(base_dir, "CATFISH_permutation_B100000_mvn_GPD", "omni_acat_mvn.csv"),
    file.path(base_dir, "catfish_omnibus_results_Fly_female_B10000",  "omni_acat_mvn.csv"),
    "/Users/nirwantandukar/Documents/Github/MAGCAT/catfish_omnibus_results_Fly_female_B10000/omni_acat_mvn.csv"
  )
  hit <- csvs[file.exists(csvs)]
  if (!length(hit)) stop("No completed CATFISH result CSV found under: ", base_dir, call. = FALSE)
  hit[[1]]
}

pick_completed_gene_table <- function(base_dir) {
  tsvs <- c(
    file.path(base_dir, "CATFISH_permutation_B100000_mvn_GPD",
              paste0(TRAIT, "_combined_genes.tsv")),
    file.path(base_dir, "CATFISH_permutation_B100000_mvn_GPD", "combined_genes.tsv")
  )
  hit <- tsvs[file.exists(tsvs)]
  if (!length(hit)) stop("No combined MAGMA gene table found under: ", base_dir, call. = FALSE)
  hit[[1]]
}

read_gene_loc_clean <- function(path) {
  x <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  req <- c("GENE", "CHR", "START", "STOP")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("Gene location file missing columns: ", paste(miss, collapse = ", "), call. = FALSE)
  x <- x[, req, drop = FALSE]
  x$GENE  <- as.character(x$GENE)
  x$CHR   <- as.integer(x$CHR)
  x$START <- as.integer(x$START)
  x$STOP  <- as.integer(x$STOP)
  x <- x[!duplicated(x$GENE), , drop = FALSE]
  x$START_EXT <- pmax(0L, x$START - WINDOW_SIZE)
  x$STOP_EXT  <- x$STOP + WINDOW_SIZE
  x
}

read_fly_gwas <- function(gwas_file) {
  x <- fread(gwas_file, header = TRUE)

  # Standardise column names
  col_map <- c("P-Value"="P","P_Value"="P","pvalue"="P","p.value"="P",
               "Positions"="POS","Position"="POS","BP"="POS","pos"="POS",
               "Chromosome"="CHR","chromosome"="CHR","chr"="CHR",
               "SNP"="SNP_ID","snp"="SNP_ID","rs"="SNP_ID")
  for (old in names(col_map)) {
    if (old %in% names(x) && !col_map[[old]] %in% names(x))
      setnames(x, old, col_map[[old]])
  }

  req <- c("CHR", "POS", "P")
  miss <- setdiff(req, names(x))
  if (length(miss)) stop("GWAS file missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

  # Map chromosome arms to numeric
  x[, CHR := CHR_MAP[as.character(CHR)]]
  x <- x[!is.na(CHR)]
  x[, CHR := as.integer(CHR)]
  x[, POS := as.integer(POS)]
  x
}

## -----------------------------------------------------------------------------
## Inputs
## -----------------------------------------------------------------------------

CATFISH_CSV     <- pick_completed_catfish_csv(BASE_DIR)
MAGMA_GENE_FILE <- pick_completed_gene_table(BASE_DIR)

cat("Using CATFISH CSV:      ", CATFISH_CSV, "\n", sep = "")
cat("Using MAGMA gene table: ", MAGMA_GENE_FILE, "\n", sep = "")
cat("GWAS file:              ", GWAS_FILE, "\n", sep = "")
cat("Output dir:             ", OUT_DIR, "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 1. Gene annotation
## -----------------------------------------------------------------------------

gene_loc <- read_gene_loc_clean(GENE_LOC_FILE)
cat("Genes in annotation after dedup: ", nrow(gene_loc), "\n", sep = "")

## -----------------------------------------------------------------------------
## 2. GWAS layer
## -----------------------------------------------------------------------------

gwas_dt <- read_fly_gwas(GWAS_FILE)
cat("Total GWAS SNPs: ", nrow(gwas_dt), "\n", sep = "")

gene_dt <- as.data.table(gene_loc)
setkey(gene_dt, CHR, START_EXT, STOP_EXT)

gwas_dt[, c("START_EXT", "STOP_EXT") := .(POS, POS)]
setkey(gwas_dt, CHR, START_EXT, STOP_EXT)

snp_gene_map <- foverlaps(
  gwas_dt, gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within", nomatch = NULL
)

gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(gwas_min_p = min(P, na.rm = TRUE), gwas_n_snps = n(), .groups = "drop") %>%
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
if (!"ZSTAT" %in% names(magma_gene)) magma_gene$ZSTAT <- NA_real_

magma_gene <- magma_gene %>%
  rename(magma_p = P, magma_z = ZSTAT) %>%
  mutate(magma_fdr = p.adjust(magma_p, method = "BH"),
         magma_bonf = p.adjust(magma_p, method = "bonferroni")) %>%
  arrange(magma_p, GENE) %>%
  mutate(magma_rank = dplyr::row_number())

if (MAGMA_MODE == "top_pct") {
  magma_n_top <- ceiling(nrow(magma_gene) * MAGMA_TOP_PCT / 100)
  magma_eff_threshold <- sort(magma_gene$magma_p)[magma_n_top]
} else {
  magma_eff_threshold <- NULL
}

cat("Genes with MAGMA results: ", nrow(magma_gene), "\n", sep = "")
if (MAGMA_MODE == "top_pct")
  cat("MAGMA effective threshold: ", format(magma_eff_threshold, digits = 4), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 4. Pathway layer
## -----------------------------------------------------------------------------

omni_results <- read.csv(CATFISH_CSV, stringsAsFactors = FALSE)

omni_col <- if ("omni_p_final" %in% names(omni_results)) "omni_p_final" else
            if ("omni_p_mvn"   %in% names(omni_results)) "omni_p_mvn"   else
            "omni_p_analytic"

# Ensure pathway_name exists
if (!"pathway_name" %in% names(omni_results)) omni_results$pathway_name <- omni_results$pathway_id

omni_results <- omni_results %>%
  arrange(.data[[omni_col]], pathway_name, pathway_id) %>%
  mutate(pathway_rank = dplyr::row_number())

pathway_genes_list <- strsplit(omni_results$genes_used, ";")
pathway_genes <- data.frame(
  pathway_id   = rep(omni_results$pathway_id,   lengths(pathway_genes_list)),
  pathway_name = rep(omni_results$pathway_name, lengths(pathway_genes_list)),
  pathway_p    = rep(omni_results[[omni_col]],   lengths(pathway_genes_list)),
  pathway_rank = rep(omni_results$pathway_rank,  lengths(pathway_genes_list)),
  GENE         = trimws(unlist(pathway_genes_list)),
  stringsAsFactors = FALSE
) %>% filter(!is.na(GENE), nzchar(GENE))

pathway_best <- pathway_genes %>%
  arrange(pathway_rank, pathway_name, pathway_id) %>%
  group_by(GENE) %>% slice_head(n = 1) %>% ungroup() %>%
  transmute(GENE, best_pathway_id = pathway_id, best_pathway_name = pathway_name,
            best_pathway_p = pathway_p, pathway_rank = pathway_rank)

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(n_top_pathways       = n_distinct(pathway_id),
            mean_pathway_mlog10p = mean(-log10(pathway_p), na.rm = TRUE),
            pathways             = paste(pathway_id,   collapse = "; "),
            pathway_names        = paste(pathway_name, collapse = "; "),
            .groups = "drop") %>%
  left_join(pathway_best, by = "GENE") %>%
  arrange(pathway_rank, GENE)

n_pathways_total <- nrow(omni_results)
cat("Total pathways used for scoring: ", n_pathways_total, "\n", sep = "")
cat("Genes with pathway membership:   ", nrow(pathway_gene_support), "\n\n", sep = "")

## -----------------------------------------------------------------------------
## 5. Integrate three layers and score genes
## -----------------------------------------------------------------------------

n_genes_total <- nrow(magma_gene)

gene_evidence <- magma_gene %>%
  left_join(gwas_gene,           by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas   = !is.na(gwas_min_p)    & gwas_min_p    <= gwas_eff_threshold,
    hit_magma  = if (MAGMA_MODE == "top_pct") {
      !is.na(magma_p) & magma_p <= magma_eff_threshold
    } else {
      !is.na(magma_fdr) & magma_fdr < MAGMA_FDR_THRESHOLD
    },
    hit_pathway    = !is.na(n_top_pathways) & n_top_pathways >= 1,
    support_layers = hit_gwas + hit_magma + hit_pathway,
    gwas_rank_score    = ifelse(!is.na(gwas_rank),    -log10(gwas_rank    / n_genes_total),    0),
    magma_rank_score   = ifelse(!is.na(magma_rank),   -log10(magma_rank   / n_genes_total),    0),
    pathway_rank_score = ifelse(!is.na(pathway_rank), -log10(pathway_rank / n_pathways_total), 0),
    multi_pathway_bonus = ifelse(!is.na(n_top_pathways), log10(n_top_pathways + 1), 0),
    score = gwas_rank_score + magma_rank_score +
            0.5 * pathway_rank_score +
            0.3 * multi_pathway_bonus
  ) %>%
  arrange(desc(score), gwas_rank, magma_rank, pathway_rank, GENE)

## -----------------------------------------------------------------------------
## 6. Save
## -----------------------------------------------------------------------------

select_cols <- c(
  "GENE", "score", "support_layers",
  "gwas_rank_score", "magma_rank_score", "pathway_rank_score", "multi_pathway_bonus",
  "magma_p", "magma_fdr", "magma_rank", "magma_z",
  "gwas_min_p", "gwas_n_snps", "gwas_rank",
  "n_top_pathways", "best_pathway_p", "pathway_rank",
  "best_pathway_id", "best_pathway_name", "pathways", "pathway_names",
  "hit_gwas", "hit_magma", "hit_pathway"
)
select_cols <- intersect(select_cols, names(gene_evidence))

top50  <- gene_evidence %>% select(all_of(select_cols)) %>% slice_head(n = 50)
top200 <- gene_evidence %>% select(all_of(select_cols)) %>% slice_head(n = 200)

all_file   <- file.path(OUT_DIR, paste0("candidate_genes_all_",   TRAIT, "_B100000_GPD.csv"))
top50_file <- file.path(OUT_DIR, paste0("candidate_genes_top50_",  TRAIT, "_B100000_GPD.csv"))
top200_file <- file.path(OUT_DIR, paste0("candidate_genes_top200_", TRAIT, "_B100000_GPD.csv"))

write.csv(gene_evidence %>% select(all_of(select_cols)), all_file,    row.names = FALSE)
write.csv(top50,                                          top50_file,  row.names = FALSE)
write.csv(top200,                                         top200_file, row.names = FALSE)

cat("Saved files:\n")
cat("  ", all_file,    "\n", sep = "")
cat("  ", top50_file,  "\n", sep = "")
cat("  ", top200_file, "\n\n", sep = "")

cat("Top 10 genes by score:\n")
print(top50 %>% slice_head(n = 10))
