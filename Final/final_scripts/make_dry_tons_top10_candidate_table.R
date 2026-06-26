#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readr)
})

window_size <- 25000L

gwas_dir <- "/Users/nirwantandukar/Documents/Research/results/GWAS/MLM/BAP/Dry_tons_per_acre"
base_dir <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre"
gene_loc_file <- "/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/sorghum.genes.loc"
out_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_supp_tables"

magma_file <- file.path(
  base_dir,
  "CATFISH_permutation_B10000_mvn_GPD_paper_tau_false",
  "Dry_tons_per_acre_combined_genes.tsv"
)

candidate_file <- file.path(
  base_dir,
  "candidate_gene_scoring_B10000_GPD_paper_tau_false",
  "candidate_genes_all_Dry_tons_per_acre_B10000_GPD_paper_tau_false.csv"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_gene_loc_clean <- function(path, window_size) {
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
  x$START_EXT <- pmax(0L, x$START - window_size)
  x$STOP_EXT <- x$STOP + window_size
  x
}

read_all_gwas <- function(dir_path) {
  files <- list.files(dir_path, pattern = "\\.assoc\\.txt$", full.names = TRUE)
  if (!length(files)) {
    stop("No GWAS .assoc.txt files found in: ", dir_path, call. = FALSE)
  }

  chr_num <- as.integer(sub(".*Chr([0-9]+)\\.maf01\\.assoc\\.txt$", "\\1", basename(files)))
  files <- files[order(chr_num)]

  pieces <- lapply(files, function(f) {
    x <- fread(f, header = TRUE)
    setnames(x, old = c("chr", "rs", "ps", "p_wald"), new = c("CHR", "SNP_ID", "POS", "P"), skip_absent = TRUE)
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
  out[, P := as.numeric(P)]
  out[]
}

format_p <- function(x) {
  ifelse(is.na(x), NA_character_, formatC(x, format = "e", digits = 2))
}

gene_loc <- read_gene_loc_clean(gene_loc_file, window_size)
gwas_dt <- read_all_gwas(gwas_dir)
magma_df <- read_tsv(magma_file, show_col_types = FALSE)
candidate_df <- read_csv(candidate_file, show_col_types = FALSE)

gene_dt <- as.data.table(gene_loc)
setkey(gene_dt, CHR, START_EXT, STOP_EXT)

top10_snps <- copy(gwas_dt[order(P, SNP_ID)][1:10])
top10_snps[, gwas_rank := seq_len(.N)]
top10_snps[, `:=`(START_EXT = POS, STOP_EXT = POS)]
setkey(top10_snps, CHR, START_EXT, STOP_EXT)

snp_gene_hits <- foverlaps(
  top10_snps,
  gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within",
  nomatch = NULL
)

top10_gwas <- as.data.frame(top10_snps) %>%
  select(CHR, SNP_ID, POS, P) %>%
  rename(gwas_chr = CHR, gwas_snp = SNP_ID, gwas_pos = POS, gwas_p = P) %>%
  left_join(
    snp_gene_hits %>%
      as.data.frame() %>%
      distinct(SNP_ID, GENE) %>%
      arrange(SNP_ID, GENE) %>%
      group_by(SNP_ID) %>%
      summarise(gwas_genes = paste(GENE, collapse = ", "), .groups = "drop"),
    by = c("gwas_snp" = "SNP_ID")
  ) %>%
  left_join(
    as.data.frame(top10_snps) %>%
      select(SNP_ID, gwas_rank) %>%
      distinct(),
    by = c("gwas_snp" = "SNP_ID")
  ) %>%
  arrange(gwas_rank, gwas_snp) %>%
  mutate(gwas_genes = ifelse(is.na(gwas_genes), "", gwas_genes))

top10_magma <- magma_df %>%
  mutate(magma_p = if ("P" %in% names(.)) P else P_MULTI) %>%
  arrange(magma_p, GENE) %>%
  transmute(
    magma_rank = row_number(),
    magma_gene = GENE,
    magma_p = magma_p
  ) %>%
  slice_head(n = 10)

top10_all3 <- candidate_df %>%
  filter(hit_gwas, hit_magma, hit_pathway) %>%
  arrange(desc(score), GENE) %>%
  transmute(
    composite_rank = row_number(),
    composite_gene = GENE,
    composite_score = round(score, 3),
    composite_best_pathway = best_pathway_name,
    composite_gwas_rank = gwas_rank,
    composite_magma_rank = magma_rank
  ) %>%
  slice_head(n = 10)

wide_table <- tibble(row_id = 1:10) %>%
  left_join(top10_gwas, by = c("row_id" = "gwas_rank")) %>%
  left_join(top10_magma, by = c("row_id" = "magma_rank")) %>%
  left_join(top10_all3, by = c("row_id" = "composite_rank"))

write_csv(top10_gwas, file.path(out_dir, "SuppTable_DryTons_top10_GWAS_SNPs_with_genes.csv"))
write_csv(top10_magma, file.path(out_dir, "SuppTable_DryTons_top10_MAGMA_genes_current_run.csv"))
write_csv(top10_all3, file.path(out_dir, "SuppTable_DryTons_top10_All3_candidate_genes_current_run.csv"))
write_csv(wide_table, file.path(out_dir, "SuppTable_DryTons_top10_candidate_panels_wide.csv"))

cat("Wrote:\n")
cat("- ", file.path(out_dir, "SuppTable_DryTons_top10_GWAS_SNPs_with_genes.csv"), "\n", sep = "")
cat("- ", file.path(out_dir, "SuppTable_DryTons_top10_MAGMA_genes_current_run.csv"), "\n", sep = "")
cat("- ", file.path(out_dir, "SuppTable_DryTons_top10_All3_candidate_genes_current_run.csv"), "\n", sep = "")
cat("- ", file.path(out_dir, "SuppTable_DryTons_top10_candidate_panels_wide.csv"), "\n\n", sep = "")

cat("Top 10 GWAS SNPs:\n")
print(top10_gwas %>% mutate(gwas_p = format_p(gwas_p)))
cat("\nTop 10 MAGMA genes:\n")
print(top10_magma %>% mutate(magma_p = format_p(magma_p)))
cat("\nTop 10 composite genes (all 3 layers):\n")
print(top10_all3)
