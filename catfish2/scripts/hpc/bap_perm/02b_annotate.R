#!/usr/bin/env Rscript
# =============================================================================
# 02b_annotate.R -- MAGMA SNP-to-gene annotation for ONE chromosome, built
# exactly as the real dry-tons analysis did (scripts/run_magma_single_chr_
# dry_tons_server.R, steps 1-3): SNP locations from the PLINK .bim, the gene
# locations with GENE/CHR/START/STOP kept and duplicate genes dropped, and
# MAGMA --annotate nonhuman window=25,25. The annotation depends only on the
# genotypes, so it is built once and every permutation reuses it.
#   Rscript 02b_annotate.R <chr2digit>
# Writes the file 03_magma_catfish.R reads: ANNOT_PATH_TEMPLATE with @CHR@.
# =============================================================================
getcfg <- function(k){ v<-Sys.getenv(k); if(v=="") stop(paste("config missing:",k)); v }
CHR        <- commandArgs(trailingOnly = TRUE)[[1]]
GENO_DIR   <- getcfg("GENO_DIR");  GENO_PAT <- getcfg("GENO_PREFIX_PATTERN")
GENE_LOC   <- getcfg("GENE_LOC_FILE")
WINDOW_KB  <- getcfg("WINDOW_KB")
MAGMA_BIN  <- getcfg("MAGMA")
annot_file <- gsub("@CHR@", CHR, getcfg("ANNOT_PATH_TEMPLATE"), fixed = TRUE)
annot_dir  <- dirname(annot_file)
prefix     <- sub("\\.genes\\.annot$", "", annot_file)
dir.create(annot_dir, recursive = TRUE, showWarnings = FALSE)
bfile <- file.path(GENO_DIR, gsub("@CHR@", CHR, GENO_PAT, fixed = TRUE))
for (f in c(paste0(bfile, ".bim"), GENE_LOC, MAGMA_BIN)) if (!file.exists(f)) stop("missing: ", f)

## 1. SNP locations from the .bim
bim <- utils::read.table(paste0(bfile, ".bim"), header = FALSE, stringsAsFactors = FALSE)
snp_loc <- data.frame(SNP = bim[[2]],
                      CHR = sub("^chr", "", as.character(bim[[1]]), ignore.case = TRUE),
                      BP  = bim[[4]], stringsAsFactors = FALSE)
snp_loc_file <- paste0(prefix, ".snp.loc.txt")
utils::write.table(snp_loc, snp_loc_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## 2. gene locations: GENE CHR START STOP, no duplicate genes, no header
x <- utils::read.table(GENE_LOC, header = TRUE, sep = "", stringsAsFactors = FALSE,
                       check.names = FALSE, comment.char = "", quote = "\"")
x <- x[, c("GENE", "CHR", "START", "STOP")]
x$CHR <- sub("^chr", "", as.character(x$CHR), ignore.case = TRUE)
x <- x[!is.na(x$GENE) & nzchar(x$GENE), ]
x <- x[!duplicated(x$GENE), ]
gene_loc_clean <- paste0(prefix, ".gene_loc.clean.txt")
utils::write.table(x, gene_loc_clean, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

## 3. MAGMA annotation
st <- system2(MAGMA_BIN, c("--annotate", "nonhuman", paste0("window=", WINDOW_KB),
                           "--snp-loc", snp_loc_file, "--gene-loc", gene_loc_clean,
                           "--out", prefix))
if (!file.exists(annot_file)) stop("MAGMA annotation did not create ", annot_file, " (exit ", st, ")")
n_genes <- sum(!grepl("^#", readLines(annot_file, warn = FALSE)))
cat(sprintf("[annot] chr%s: %d SNPs, %d genes annotated -> %s\n", CHR, nrow(snp_loc), n_genes, annot_file))
