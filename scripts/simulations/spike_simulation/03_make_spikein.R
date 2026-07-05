#!/usr/bin/env Rscript
# =============================================================================
# 12_make_spikein.R  --  Build staged spike-in synthetic phenotypes on the REAL
#                        genotypes using the extracted causal dosage matrix.
# =============================================================================
suppressMessages({
  library(data.table)
})

getcfg <- function(k) {
  v <- Sys.getenv(k, "")
  if (!nzchar(v)) stop(paste("config missing:", k), call. = FALSE)
  v
}

SPIKE_WORK <- getcfg("SPIKE_WORK")
PHENO_DIR <- getcfg("PHENO_DIR")
NTOTAL_FILE <- getcfg("NTOTAL_FILE")
GENO_DIR <- getcfg("GENO_DIR")
GENO_PAT <- getcfg("GENO_PREFIX_PATTERN")
TRAIT <- getcfg("TRAIT")
ARCH <- getcfg("SPIKE_ARCH")
H2 <- as.numeric(getcfg("SPIKE_H2"))
N_SIM <- as.integer(getcfg("SPIKE_N_SIM"))
SEED <- as.integer(getcfg("SPIKE_SEED"))
N_STRONG <- as.integer(getcfg("SPIKE_N_STRONG"))
MSA_FRAC <- as.numeric(getcfg("SPIKE_MSA_FRAC"))
PATHWAY <- getcfg("SPIKE_PATHWAY")

raw <- fread(file.path(SPIKE_WORK, "causal_geno.raw"))
iid <- as.character(raw$IID)
geno_cols <- setdiff(names(raw), c("FID", "IID", "PAT", "MAT", "SEX", "PHENOTYPE"))
G <- as.matrix(raw[, ..geno_cols])
for (j in seq_len(ncol(G))) {
  cj <- G[, j]
  m <- mean(cj, na.rm = TRUE)
  cj[is.na(cj)] <- m
  s <- sd(cj)
  if (!is.finite(s) || s == 0) s <- 1
  G[, j] <- (cj - m) / s
}

cmap <- fread(file.path(SPIKE_WORK, "causal_map.tsv"))
strip <- function(x) sub("_[ACGT]+$", "", x)
col_snp <- strip(geno_cols)
gene_of_col <- cmap$GENE[match(col_snp, cmap$SNP)]
ok <- !is.na(gene_of_col)
G <- G[, ok, drop = FALSE]
gene_of_col <- gene_of_col[ok]
genes <- unique(gene_of_col)
K <- length(genes)
cat(sprintf("[spike] %d causal genes usable, %d individuals\n", K, nrow(G)))

beta_gene <- setNames(numeric(K), genes)
set.seed(SEED)
if (ARCH == "SGP") {
  beta_gene[sample(genes, 1)] <- 1
} else if (ARCH == "SSS") {
  beta_gene[sample(genes, min(N_STRONG, K))] <- 1
} else if (ARCH == "MDS" || ARCH == "WDE") {
  beta_gene[] <- 1
} else if (ARCH == "MSA") {
  big <- sample(genes, 1)
  beta_gene[] <- sqrt((1 - MSA_FRAC) / max(K - 1, 1))
  beta_gene[big] <- sqrt(MSA_FRAC)
} else {
  stop(paste("unknown SPIKE_ARCH:", ARCH), call. = FALSE)
}
beta_col <- beta_gene[gene_of_col]

fam_pfx <- file.path(GENO_DIR, gsub("@CHR@", "01", GENO_PAT, fixed = TRUE))
fam <- fread(paste0(fam_pfx, ".fam"), header = FALSE)
fam_iids <- as.character(fam$V2)
row_idx <- match(fam_iids, iid)

dir.create(PHENO_DIR, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(sum(!is.na(row_idx))), NTOTAL_FILE)

manifest <- data.table(perm = integer(), seed = integer(), file = character())
for (b in seq_len(N_SIM)) {
  seed_b <- SEED + b
  set.seed(seed_b)
  V <- as.numeric(G %*% beta_col)
  V <- (V - mean(V)) / sd(V)
  e <- rnorm(length(V))
  y <- sqrt(H2) * V + sqrt(1 - H2) * e
  y_fam <- y[row_idx]
  out <- data.table(ID = fam_iids, TRAIT = y_fam)
  setnames(out, c("ID", TRAIT))
  outf <- file.path(PHENO_DIR, sprintf("pheno_perm_%03d.csv", b))
  fwrite(out, outf, na = "NA")
  manifest <- rbind(manifest, data.table(perm = b, seed = seed_b, file = outf))
}
fwrite(manifest, file.path(PHENO_DIR, "manifest.tsv"), sep = "\t")
fwrite(
  data.table(
    spiked_pathway = PATHWAY,
    archetype = ARCH,
    h2 = H2,
    n_causal_genes = K,
    causal_genes = paste(genes, collapse = ";")
  ),
  file.path(PHENO_DIR, "truth.tsv"),
  sep = "\t"
)
cat(sprintf("[spike] wrote %d synthetic phenotype CSVs (arch=%s, h2=%.3f) + truth.tsv\n",
            N_SIM, ARCH, H2))
