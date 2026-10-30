#!/usr/bin/env Rscript
## =============================================================================
## 20_lipid_manhattan_prep.R -- thin the SNP-level GLGC results for Manhattan
## plots (Figure 13, drawn by 21_lipid_results.R).
##
## The unfloored -log10(p) files carry SNP and NEGLOG10P only, ~20 million rows
## per trait. Positions come from the reference-panel .bim (the same panel MAGMA
## used). Every SNP at or above KEEP_ABOVE is kept; below it a random fraction
## FRAC_BELOW is kept, which preserves the look of the plot at ~3% of the size.
##
## Usage (from catfish2/scripts):
##   Rscript 20_lipid_manhattan_prep.R
##   TRAITS=hdl Rscript 20_lipid_manhattan_prep.R          # one trait only
##   GWAS_FILE_HDL=/path/to/hdl_eur.neglog10.txt Rscript 20_lipid_manhattan_prep.R
##
## Writes ../results/20_manhattan_snps_<trait>.tsv  (CHR BP SNP NEGLOG10P;
## uncompressed, because data.table on some installs has no zlib for fwrite)
## Needs ~4 GB of memory while the .bim and one trait are loaded.
## =============================================================================
suppressPackageStartupMessages(library(data.table))

DATA_D     <- Sys.getenv("CATFISH_DATA", path.expand("~/Documents/Github/MAGCAT_data"))
BIM        <- Sys.getenv("BIM", file.path(DATA_D, "g1000_eur.bim"))
TRAITS     <- strsplit(Sys.getenv("TRAITS", "ldl,hdl"), ",")[[1]]
KEEP_ABOVE <- as.numeric(Sys.getenv("KEEP_ABOVE", "3"))
FRAC_BELOW <- as.numeric(Sys.getenv("FRAC_BELOW", "0.01"))
OUT_DIR    <- Sys.getenv("OUT_DIR", "../results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

find_gwas <- function(t) {
  ov <- Sys.getenv(sprintf("GWAS_FILE_%s", toupper(t)), "")
  if (nzchar(ov)) return(ov)
  cand <- file.path(DATA_D, c(sprintf("lipids/%s_eur.neglog10.txt", t),
                              sprintf("%s/%s_eur.neglog10.txt", t, t),
                              sprintf("%s_eur.neglog10.txt", t)))
  hit <- cand[file.exists(cand)]
  if (!length(hit))
    stop(sprintf("no %s_eur.neglog10.txt under %s (looked in lipids/, %s/ and the top level); set GWAS_FILE_%s",
                 t, DATA_D, t, toupper(t)), call. = FALSE)
  hit[1]
}

if (!file.exists(BIM)) stop("reference .bim not found: ", BIM, call. = FALSE)
cat("reading", BIM, "\n")
bim <- fread(BIM, header = FALSE, select = c(1, 2, 4),
             col.names = c("CHR", "SNP", "BP"), showProgress = FALSE)
bim <- bim[CHR %in% as.character(1:22) | CHR %in% 1:22]
bim[, CHR := as.integer(CHR)]
setkey(bim, SNP)
cat(sprintf("  %d autosomal reference SNPs\n", nrow(bim)))

for (t in TRAITS) {
  f <- find_gwas(t)
  cat(sprintf("[%s] reading %s\n", t, f))
  nl <- fread(f, select = c("SNP", "NEGLOG10P"), showProgress = FALSE)
  m <- bim[nl, on = "SNP", nomatch = 0L]
  m <- m[is.finite(NEGLOG10P)]
  set.seed(20261029)
  keep <- m$NEGLOG10P >= KEEP_ABOVE | stats::runif(nrow(m)) < FRAC_BELOW
  m <- m[keep][order(CHR, BP)]
  out <- file.path(OUT_DIR, sprintf("20_manhattan_snps_%s.tsv", t))
  fwrite(m[, list(CHR, BP, SNP, NEGLOG10P)], out, sep = "\t")
  cat(sprintf("  %d SNPs with position; kept %d (all with -log10 p >= %g); max -log10 p %.1f\n  wrote %s\n",
              nrow(nl), nrow(m), KEEP_ABOVE, max(m$NEGLOG10P), out))
  rm(nl, m); invisible(gc())
}
