#!/usr/bin/env Rscript

suppressMessages({ library(data.table) })
data.table::setDTthreads(1L)

getcfg <- function(k){ v <- Sys.getenv(k); if (v == "") stop(paste("config missing:", k)); v }

N_PERM <- as.integer(getcfg("N_PERM"))
MASTER_SEED <- as.integer(getcfg("MASTER_SEED"))
TRAIT <- getcfg("TRAIT")
PHENO_CSV <- getcfg("PHENO_CSV")
PHENO_DIR <- getcfg("PHENO_DIR")
NTOTAL_FILE <- getcfg("NTOTAL_FILE")
GENO_DIR <- getcfg("GENO_DIR")
GENO_PAT <- getcfg("GENO_PREFIX_PATTERN")

fam_pfx <- file.path(GENO_DIR, gsub("@CHR@", "01", GENO_PAT, fixed = TRUE))
fam_file <- paste0(fam_pfx, ".fam")
if (!file.exists(fam_file)) {
  stop(
    "Missing PLINK family file: ", fam_file, "\n",
    "Run 00_prepare_genotype.sh first and confirm VCF_DIR / VCF_PREFIX_PATTERN ",
    "point to the real per-chromosome SAP VCF files.",
    call. = FALSE
  )
}
fam <- fread(fam_file, header = FALSE)
geno_iids <- as.character(fam$V2)

ph <- fread(PHENO_CSV)
id_col <- names(ph)[1]
if (!TRAIT %in% names(ph)) stop(sprintf("Trait '%s' not in %s", TRAIT, PHENO_CSV))
ph[[id_col]] <- as.character(ph[[id_col]])

val <- ph[[TRAIT]][match(geno_iids, ph[[id_col]])]
n_obs <- sum(!is.na(val))
cat(sprintf("[perm] genotyped individuals: %d | non-missing %s: %d\n",
            length(geno_iids), TRAIT, n_obs))
if (n_obs < 50) stop("Almost nothing matched between CSV IDs and .fam IDs.")

dir.create(PHENO_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(n_obs), NTOTAL_FILE)

obs <- which(!is.na(val))
manifest <- data.table(perm = integer(), seed = integer(), file = character())
for (b in seq_len(N_PERM)) {
  seed_b <- MASTER_SEED + b
  set.seed(seed_b)
  pv <- val
  pv[obs] <- sample(val[obs])
  out <- data.table(ID = geno_iids, TRAIT = pv)
  setnames(out, c(id_col, TRAIT))
  outf <- file.path(PHENO_DIR, sprintf("pheno_perm_%03d.csv", b))
  fwrite(out, outf, na = "NA")
  manifest <- rbind(manifest, data.table(perm = b, seed = seed_b, file = outf))
}
fwrite(manifest, file.path(PHENO_DIR, "manifest.tsv"), sep = "\t")
cat(sprintf("[perm] wrote %d permuted phenotype CSVs to %s\n", N_PERM, PHENO_DIR))
