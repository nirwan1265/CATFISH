#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(CATFISH)
})
data.table::setDTthreads(1L)

getcfg <- function(k) {
  v <- Sys.getenv(k)
  if (!nzchar(v)) stop("config missing: ", k, call. = FALSE)
  v
}

options(magma.path = getcfg("MAGMA"))

GENE_LOC_FILE <- getcfg("GENE_LOC_FILE")
GENO_DIR <- getcfg("GENO_DIR")
ANNOT_DIR <- getcfg("ANNOT_DIR")
GENO_PREFIX_PATTERN <- getcfg("GENO_PREFIX_PATTERN")
CHROMS <- strsplit(trimws(getcfg("CHROMS")), "\\s+")[[1]]
PLINK <- getcfg("PLINK")
VCF_DIR <- getcfg("VCF_DIR")
VCF_PREFIX_PATTERN <- getcfg("VCF_PREFIX_PATTERN")

dir.create(GENO_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(ANNOT_DIR, recursive = TRUE, showWarnings = FALSE)

for (chr in CHROMS) {
  chr_num <- sub("^0+", "", chr)
  vcf_file <- file.path(
    VCF_DIR,
    gsub("@CHR_NUM@", chr_num, gsub("@CHR@", chr, VCF_PREFIX_PATTERN, fixed = TRUE), fixed = TRUE)
  )
  if (!file.exists(vcf_file)) stop("missing VCF file: ", vcf_file, call. = FALSE)

  plink_prefix <- file.path(GENO_DIR, gsub("@CHR@", chr, GENO_PREFIX_PATTERN, fixed = TRUE))
  if (!file.exists(paste0(plink_prefix, ".bed"))) {
    cat("[prep] VCF -> PLINK for chr", chr, "\n")
    status <- system2(PLINK, c(
      "--vcf", vcf_file,
      "--make-bed",
      "--out", plink_prefix,
      "--allow-no-sex",
      "--allow-extra-chr"
    ))
    if (!identical(status, 0L)) stop("plink --make-bed failed for chr ", chr, call. = FALSE)
  }
}

for (chr_label in CHROMS) {
  chr_num <- as.integer(sub("^0+", "", chr_label))
  chr_out_dir <- file.path(ANNOT_DIR, paste0("chr", chr_label))
  dir.create(chr_out_dir, recursive = TRUE, showWarnings = FALSE)
  out_prefix <- sprintf("stem_volume_chr%s", chr_label)
  annot_file <- file.path(chr_out_dir, paste0(out_prefix, ".genes.annot"))
  if (file.exists(annot_file)) next

  plink_prefix <- file.path(GENO_DIR, gsub("@CHR@", chr_label, GENO_PREFIX_PATTERN, fixed = TRUE))
  bim_file <- paste0(plink_prefix, ".bim")
  if (!file.exists(bim_file)) {
    stop("Missing PLINK BIM file for chr", chr_label, ": ", bim_file, call. = FALSE)
  }
  bim <- fread(bim_file, header = FALSE)
  if (ncol(bim) < 4L) stop("Unexpected BIM format in: ", bim_file, call. = FALSE)
  snp_stats <- data.table(
    chr = as.integer(bim[[1]]),
    rs = as.character(bim[[2]]),
    ps = as.integer(bim[[4]]),
    p_wald = 1
  )
  if (!nrow(snp_stats)) {
    stop("No SNPs found in BIM for chr", chr_label, call. = FALSE)
  }
  tmp_stats <- file.path(chr_out_dir, sprintf("stem_volume_chr%s_stats.txt", chr_label))
  fwrite(snp_stats, tmp_stats, sep = "\t")
  ann <- magma_annotate(
    stats_file = tmp_stats,
    rename_columns = c(CHR = "chr", SNP = "rs", POS = "ps", PVALUE = "p_wald"),
    gene_loc = GENE_LOC_FILE,
    out_prefix = out_prefix,
    out_dir = chr_out_dir,
    window = c(10, 10),
    sep = "\t",
    nonhuman = TRUE
  )
  cat("[prep] chr", chr_label, " annotate target:", ann$gene_annot, "\n")
  if (!file.exists(annot_file)) {
    stop(
      "MAGMA annotation failed to create: ", annot_file, "\n",
      "Check MAGMA stderr/stdout from 00_prepare_genotype.sh and verify GENE_LOC_FILE/MAGMA path.",
      call. = FALSE
    )
  }
}

cat("[prep] genotype and annotation preparation complete\n")
