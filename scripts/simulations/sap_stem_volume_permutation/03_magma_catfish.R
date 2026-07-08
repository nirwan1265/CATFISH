#!/usr/bin/env Rscript

suppressPackageStartupMessages({ library(CATFISH) })
getcfg <- function(k){ v <- Sys.getenv(k); if (v == "") stop(paste("config missing:",k)); v }

args <- commandArgs(trailingOnly = TRUE)
B <- as.integer(args[[1]])
BPAD <- sprintf("%03d", B)
CHR <- args[[2]]

GENO_DIR <- getcfg("GENO_DIR")
GENO_PAT <- getcfg("GENO_PREFIX_PATTERN")
ANNOT_TMPL <- getcfg("ANNOT_PATH_TEMPLATE")
GWAS_DIR <- getcfg("GWAS_DIR")
MAGMA_OUT <- getcfg("MAGMA_OUT")
MAGMA_BIN <- getcfg("MAGMA")
GENE_MODEL <- getcfg("GENE_MODEL")
NTOTAL_FILE <- getcfg("NTOTAL_FILE")

N_TOTAL <- as.integer(readLines(NTOTAL_FILE, warn = FALSE)[1])
if (is.na(N_TOTAL) || N_TOTAL <= 0) stop("bad N_TOTAL from ", NTOTAL_FILE)

BFILE_PREFIX <- file.path(GENO_DIR, gsub("@CHR@", CHR, GENO_PAT, fixed = TRUE))
gene_annot_file <- gsub("@CHR@", CHR, ANNOT_TMPL, fixed = TRUE)
GWAS_FILE <- file.path(GWAS_DIR, sprintf("perm_%s", BPAD), sprintf("chr%s", CHR), "assoc.txt")
MAGMA_DIR <- file.path(MAGMA_OUT, sprintf("perm_%s", BPAD), sprintf("chr%s", CHR))
OUT_PREFIX <- sprintf("perm_%s_chr%s", BPAD, CHR)
dir.create(MAGMA_DIR, recursive = TRUE, showWarnings = FALSE)

for (f in c(GWAS_FILE, gene_annot_file, paste0(BFILE_PREFIX, ".bed"))) {
  if (!file.exists(f)) stop("missing input: ", f)
}

find_first_match <- function(choices, available) {
  hit <- choices[choices %in% available]
  if (length(hit)) hit[[1]] else NA_character_
}

hdr <- utils::read.table(
  GWAS_FILE, header = TRUE, sep = "", nrows = 5,
  stringsAsFactors = FALSE, check.names = FALSE,
  comment.char = "", quote = "\""
)
cn <- colnames(hdr)
RENAME_COLUMNS <- c(
  SNP = find_first_match(c("SNP","rs","RS","Marker","marker","ID","id"), cn),
  PVALUE = find_first_match(c("P","P.value","P_VALUE","p_wald","P_WALD","pvalue","PVAL"), cn),
  NMISS = find_first_match(c("NMISS","n_miss","nmiss","N_MISS"), cn)
)
if (any(is.na(RENAME_COLUMNS))) {
  stop("could not detect SNP/P/NMISS in ", GWAS_FILE, " | cols: ", paste(cn, collapse = ", "))
}

options(magma.path = MAGMA_BIN)
CATFISH::magma_gene(
  bfile = BFILE_PREFIX,
  gene_annot = gene_annot_file,
  stats_file = GWAS_FILE,
  sep = "",
  n_total = N_TOTAL,
  rename_columns = RENAME_COLUMNS,
  out_prefix = OUT_PREFIX,
  out_dir = MAGMA_DIR,
  gene_model = GENE_MODEL,
  n_threads = 1
)

model_tag <- gsub("[^A-Za-z0-9]+", "_", GENE_MODEL)
cat(sprintf("[magma] wrote %s\n",
            file.path(MAGMA_DIR, paste0(OUT_PREFIX, ".", model_tag, ".genes.out"))))
