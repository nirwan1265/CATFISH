#!/usr/bin/env Rscript
## =============================================================================
## 02_simulate_z.R <cell> -- standardise plink2 --variant-score output to z.
## plink2 gives S_j = sum_i g_ij u_i (g = ALT dosage). With m_j and s_j the mean
## and SD of g_j over the panel's people (from --geno-counts, step 00),
##     z_j = (S_j - m_j sum_i u_i) / s_j,
## which is X'u for standardised genotypes. Checked against step 01's exact z at
## the causal SNPs; the job stops if they disagree.
## Then adds the signal: each causal SNP's lambda spread to the SNPs within
## LD_WINDOW_KB through their panel correlation (local LD only; see step 01).
## Writes RUN/z/z_cell<k>.tsv.gz (SNP, CHR, rep_...) and removes the vscore file.
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
K <- commandArgs(trailingOnly = TRUE)[1]; RUN <- env("RUN")
vs_file <- list.files(file.path(RUN, "z"), pattern = sprintf("^vscore_cell%s\\.vscore", K), full.names = TRUE)
if (length(vs_file) != 1) stop("expected one plink2 .vscore file for cell ", K)
vs <- fread(vs_file)
setnames(vs, sub("^#", "", names(vs)))
reps <- grep("^rep_", names(vs), value = TRUE)
if (!length(reps) || !"ID" %in% names(vs)) stop("unexpected .vscore columns: ", paste(names(vs), collapse = " "))

gc <- fread(env("SIM_COUNTS")); setnames(gc, sub("^#", "", names(gc)))
n <- gc$HOM_REF_CT + gc$HET_REF_ALT_CTS + gc$TWO_ALT_GENO_CTS
m <- (gc$HET_REF_ALT_CTS + 2 * gc$TWO_ALT_GENO_CTS) / n
s <- sqrt((gc$HET_REF_ALT_CTS + 4 * gc$TWO_ALT_GENO_CTS) / n - m^2)
i <- match(vs$ID, gc$ID); if (anyNA(i)) stop("vscore SNPs missing from the genotype counts")

U <- fread(file.path(RUN, "u", sprintf("u_cell%s.tsv", K)))
sum_u <- colSums(as.matrix(U[, ..reps]))
Z <- (as.matrix(vs[, ..reps]) - outer(m[i], sum_u)) / s[i]

## check against the exact z of the causal SNPs (none in S0)
zc <- fread(file.path(RUN, "u", sprintf("zcheck_cell%s.tsv", K)))
if (nrow(zc)) {
  zc[, z_plink := Z[cbind(match(snp, vs$ID), match(sprintf("rep_%04d", rep_id), reps))]]
  ## step 01 reads genotypes with plink2 --export A, step 02 with --variant-score
  ## (ALT dosage). If the two count different alleles, a SNP's z changes sign
  ## only, which does not matter here (causal effects have random sign); a
  ## wrong sample order or SNP mix-up changes the values and stops the job.
  zc[, d_same := abs(z_plink - z_exact)][, d_flip := abs(z_plink + z_exact)]
  zc[, d := pmin(d_same, d_flip)]
  fwrite(zc, file.path(RUN, "z", sprintf("zcheck_result_cell%s.tsv", K)), sep = "\t")
  cat(sprintf("[z] cell %s: %d causal SNPs; max |z_plink - z_exact| = %.2e, max |z_plink + z_exact| = %.2e, %d of them sign-flipped; max difference up to sign = %.2e\n",
              K, nrow(zc), max(zc$d_same), max(zc$d_flip), sum(zc$d_flip < zc$d_same), max(zc$d)))
  print(head(zc[order(-d)], 5))
  if (!is.finite(max(zc$d)) || max(zc$d) > 1e-3)
    stop("plink2 z does not match the exact z beyond a sign (sample order or SNP IDs); see ",
         file.path(RUN, "z", sprintf("zcheck_result_cell%s.tsv", K)))
}
cat(sprintf("[z] cell %s: noise only: SD of z per replicate %.3f-%.3f (about 1)\n",
            K, min(apply(Z, 2, sd)), max(apply(Z, 2, sd))))

## ---- signal: each causal SNP's lambda spread through local LD ------------------
## mean z_j = sum_k r_jk lambda_k over causal SNPs k within LD_WINDOW_KB of j, with
## r the panel correlation (both SNPs read with the same allele coding, so r does
## not depend on which allele is counted; lambda has random sign)
cs <- fread(file.path(RUN, "causal.tsv"))
cs <- cs[rep_id %in% as.integer(sub("rep_", "", reps))]
if (nrow(cs)) {
  WIN <- as.numeric(env("LD_WINDOW_KB")) * 1000; PANEL <- env("SIM_PANEL"); PLINK2 <- env("PLINK2")
  bim <- fread(paste0(PANEL, ".bim"), header = FALSE, col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
  row_of <- match(bim$snp, vs$ID)                  # bim order -> Z rows
  for (cc in sort(unique(cs$chr))) {
    pre <- tempfile()
    st <- system2(PLINK2, c("--bfile", PANEL, "--chr", cc, "--export", "A", "--out", pre),
                  stdout = FALSE, stderr = FALSE)
    if (st != 0) stop("plink2 --export A failed for chromosome ", cc)
    G <- as.matrix(fread(paste0(pre, ".raw"), drop = 1:6)); unlink(paste0(pre, c(".raw", ".log")))
    bi <- which(bim$chr == cc)
    if (ncol(G) != length(bi)) stop("chromosome ", cc, ": exported SNPs do not match the panel")
    m <- colMeans(G); sdv <- sqrt(colMeans(G^2) - m^2)
    X <- sweep(sweep(G, 2, m, "-"), 2, sdv, "/"); rm(G)
    pos <- bim$pos[bi]; n <- nrow(X)
    for (q in which(cs$chr == cc)) {
      k <- match(cs$snp[q], bim$snp[bi])
      w <- which(abs(pos - pos[k]) <= WIN)
      r <- as.vector(crossprod(X[, w, drop = FALSE], X[, k])) / n
      col <- match(sprintf("rep_%04d", cs$rep_id[q]), reps)
      rows <- row_of[bi[w]]
      Z[rows, col] <- Z[rows, col] + cs$lambda[q] * r
    }
    rm(X); gc(verbose = FALSE)
  }
  cat(sprintf("[z] cell %s: signal from %d causal SNPs added through LD within %g kb\n",
              K, nrow(cs), WIN / 1000))
}
cat(sprintf("[z] cell %s: %d SNPs x %d replicates; SD of z per replicate %.3f-%.3f\n",
            K, nrow(Z), ncol(Z), min(apply(Z, 2, sd)), max(apply(Z, 2, sd))))
out <- cbind(data.table(SNP = vs$ID, CHR = gc$CHROM[i]), as.data.table(round(Z, 6)))
fwrite(out, file.path(RUN, "z", sprintf("z_cell%s.tsv.gz", K)), sep = "\t")
file.remove(vs_file)
cat(sprintf("[z] wrote %s\n", file.path(RUN, "z", sprintf("z_cell%s.tsv.gz", K))))
