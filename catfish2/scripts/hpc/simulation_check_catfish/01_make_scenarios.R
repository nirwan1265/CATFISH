#!/usr/bin/env Rscript
## =============================================================================
## 01_make_scenarios.R -- the design, the truth, and the per-person vectors u
## from which step 02 makes each replicate's GWAS statistics.
##
## For standardised genotypes X (n = 503 people x SNPs) and u = e / sqrt(n),
## e ~ N(0, I_n), the noise z = X'u has exactly the panel's LD as covariance.
## Step 02 computes X'u for every panel SNP and adds the signal: each causal
## SNP's expected z (lambda, random sign) spread to the SNPs within
## LD_WINDOW_KB through their panel correlation. Signal is spread through local
## LD only: with 503 people, unlinked SNPs show chance correlations of about
## +/- 0.045, and spreading through those would inflate every SNP.
## Here u is drawn and the causal SNPs' noise z is computed exactly (from their
## genotypes) for step 02's check.
##
## Scenarios (archetype layers as in 17_archetype_power_wy.R: in a pathway of m
## genes, a layer marks max(1, round(m x prop)) genes, disjoint across layers,
## each causal gene getting one causal SNP with |lambda| = delta x EFFECT_SCALE):
##   S0  no causal SNPs
##   S1  background only: BG_FRAC of the genes, one causal SNP each at BG_LAMBDA
##   S2 SSS 15% x 3.00 | S3 MDS 60% x 1.10 | S4 WDE 90% x 0.55
##   S5 MSA 10% x 1.80 + 60% x 1.10 | S6 SGP 5% x 4.00   (each + background)
## Target pathways: one fixed non-lipid Reactome set per size class (about 15,
## 40 and 100 genes on chr 1-10). Background genes exclude the target's genes.
##
## MODE=pilot: S0, and S3 on the medium target at each of PILOT_SCALES.
## MODE=full : S0, S1, and S2-S6 x three sizes at EFFECT_SCALE.
## Writes RUN/design.tsv, RUN/targets.tsv, RUN/causal.tsv,
##        RUN/u/u_cell<k>.tsv (#FID IID rep_...), RUN/u/zcheck_cell<k>.tsv
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
MODE <- env("MODE"); RUN <- env("RUN"); IN <- env("IN")
SEED <- as.integer(env("MASTER_SEED")); PLINK2 <- env("PLINK2"); PANEL <- env("SIM_PANEL")
BG_FRAC <- as.numeric(env("BG_FRAC")); BG_LAMBDA <- as.numeric(env("BG_LAMBDA"))
dir.create(file.path(RUN, "u"), recursive = TRUE, showWarnings = FALSE)

ARCH <- list(
  SSS = list(c(prop = 0.15, delta = 3.00)),
  MDS = list(c(prop = 0.60, delta = 1.10)),
  WDE = list(c(prop = 0.90, delta = 0.55)),
  MSA = list(c(prop = 0.10, delta = 1.80), c(prop = 0.60, delta = 1.10)),
  SGP = list(c(prop = 0.05, delta = 4.00)))
SCODE <- c(SSS = "S2", MDS = "S3", WDE = "S4", MSA = "S5", SGP = "S6")
SIZES <- list(small = c(12, 18), medium = c(35, 45), large = c(90, 110))
LIPID <- "LIPID|LIPOPROTEIN|CHOLESTEROL|LDL|HDL|CHYLOMICRON|BILE|STEROL|STEROID|FATTY|TRIGLYCERIDE|TRIACYLGLYCEROL|LIPASE|LPL|PPAR|NR1H|SREBP|SREBF|PHOSPHOLIPID|SPHINGO|FAT_SOLUBLE|KETONE|ACYL"

## ---- genes, SNPs, pathways --------------------------------------------------------
genes <- fread(file.path(IN, "genes.tsv"), colClasses = list(character = "gene"))
gsnp <- list()
for (x in readLines(env("SIM_ANNOT"))) {
  if (startsWith(x, "#")) next
  f <- strsplit(x, "\t", fixed = TRUE)[[1]]; gsnp[[f[1]]] <- f[-(1:2)]
}
bim <- fread(paste0(PANEL, ".bim"), header = FALSE, col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
sets <- list()
for (x in readLines(env("SIM_SETS"))) { f <- strsplit(x, "\t", fixed = TRUE)[[1]]; sets[[f[1]]] <- f[-1] }
cat(sprintf("[design] %d genes, %d panel SNPs, %d pathways\n", nrow(genes), nrow(bim), length(sets)))

## ---- targets: one non-lipid pathway per size class (fixed by MASTER_SEED) ---------
set.seed(SEED)
ns <- vapply(sets, length, integer(1))
targets <- rbindlist(lapply(names(SIZES), function(sz) {
  r <- SIZES[[sz]]
  cand <- names(sets)[ns >= r[1] & ns <= r[2] & !grepl(LIPID, names(sets))]
  if (!length(cand)) stop("no non-lipid pathway with ", r[1], "-", r[2], " genes")
  id <- cand[sample.int(length(cand), 1)]
  data.table(size_class = sz, target_id = id, target_n = length(sets[[id]]))
}))
fwrite(targets, file.path(RUN, "targets.tsv"), sep = "\t")
print(targets)

## ---- design -----------------------------------------------------------------------
cells <- if (MODE == "pilot") {
  scales <- as.numeric(strsplit(env("PILOT_SCALES"), " ")[[1]])
  rbind(data.table(scenario = "S0", archetype = "none", size_class = NA_character_, scale = 0),
        data.table(scenario = "S3", archetype = "MDS", size_class = "medium", scale = scales))
} else {
  sc <- as.numeric(env("EFFECT_SCALE"))
  rbind(data.table(scenario = c("S0", "S1"), archetype = c("none", "background"),
                   size_class = NA_character_, scale = 0),
        CJ(archetype = names(ARCH), size_class = names(SIZES), sorted = FALSE)[
          , .(scenario = SCODE[archetype], archetype, size_class, scale = sc)])
}
NREP <- as.integer(env(if (MODE == "pilot") "N_REP_PILOT" else "N_REP_FULL"))
cells[, cell := .I]
design <- cells[rep(seq_len(.N), each = NREP)][, `:=`(rep_id = .I, col = rowid(cell))]
design <- merge(design, targets, by = "size_class", all.x = TRUE, sort = FALSE)[order(rep_id)]
setcolorder(design, c("rep_id", "cell", "col", "scenario", "archetype", "size_class", "scale"))
fwrite(design, file.path(RUN, "design.tsv"), sep = "\t")
cat(sprintf("[design] %s: %d cells, %d replicates\n", MODE, nrow(cells), nrow(design)))

## ---- causal SNPs per replicate ----------------------------------------------------
## a SNP ID that occurs twice in the panel cannot be matched unambiguously
## between steps, so it is never chosen as a causal SNP
dup_ids <- bim$snp[duplicated(bim$snp)]
cat(sprintf("[design] %d SNP IDs occur more than once in the panel; not used as causal SNPs\n", length(unique(dup_ids))))
pick_snp <- function(g) { s <- setdiff(gsnp[[g]], dup_ids); if (!length(s)) NA_character_ else s[sample.int(length(s), 1)] }
causal <- rbindlist(lapply(seq_len(nrow(design)), function(i) {
  d <- design[i]; set.seed(SEED + d$rep_id)
  out <- list(); tg <- character(0)
  if (d$archetype %in% names(ARCH)) {
    tg <- sets[[d$target_id]]; m <- length(tg); avail <- sample(tg)
    start <- 1L
    for (ly in ARCH[[d$archetype]]) {
      n_ly <- max(1L, as.integer(round(m * ly[["prop"]])))
      gg <- avail[start:(start + n_ly - 1L)]; start <- start + n_ly
      out[[length(out) + 1]] <- data.table(gene = gg, role = "target",
                                           abs_lambda = ly[["delta"]] * d$scale)
    }
  }
  if (d$scenario != "S0") {
    pool <- setdiff(genes$gene, tg)
    gg <- pool[sample.int(length(pool), round(BG_FRAC * nrow(genes)))]
    out[[length(out) + 1]] <- data.table(gene = gg, role = "background", abs_lambda = BG_LAMBDA)
  }
  if (!length(out)) return(NULL)
  o <- rbindlist(out)
  o[, snp := vapply(gene, pick_snp, character(1))]
  o <- o[!is.na(snp)]                             # gene with only duplicated SNP IDs: skipped
  o[, lambda := abs_lambda * sample(c(-1, 1), .N, replace = TRUE)]
  o[, rep_id := d$rep_id][]
}))
if (nrow(causal)) {
  causal <- merge(causal, bim[, .(snp, chr, pos)], by = "snp", sort = FALSE)
  ## two causal genes can draw the same SNP: their effects add
  causal <- causal[, .(gene = paste(gene, collapse = ","), role = paste(unique(role), collapse = ","),
                       lambda = sum(lambda)), by = .(rep_id, snp, chr, pos)]
}
fwrite(causal, file.path(RUN, "causal.tsv"), sep = "\t")
cat(sprintf("[design] %d causal SNP entries\n", nrow(causal)))

## ---- u for every replicate, one file per cell --------------------------------------
fam <- fread(paste0(PANEL, ".fam"), header = FALSE, select = 1:2, col.names = c("FID", "IID"))
n <- nrow(fam)
for (k in cells$cell) {
  reps <- design[cell == k, rep_id]
  cs <- causal[rep_id %in% reps]
  X <- NULL
  if (nrow(cs)) {
    ids <- unique(cs$snp)
    lst <- tempfile(fileext = ".txt"); writeLines(ids, lst); pre <- tempfile()
    st <- system2(PLINK2, c("--bfile", PANEL, "--extract", lst, "--export", "A",
                            "--out", pre), stdout = FALSE, stderr = FALSE)
    if (st != 0) stop("plink2 --export A failed for cell ", k)
    A <- fread(paste0(pre, ".raw"))
    ord <- bim$snp[bim$snp %in% ids]                 # export order = panel order
    G <- as.matrix(A[, -(1:6)]); colnames(G) <- ord
    stopifnot(nrow(G) == n, identical(as.character(A$IID), as.character(fam$IID)))
    m <- colMeans(G); s <- sqrt(colMeans(G^2) - m^2)
    if (any(!is.finite(s) | s <= 0)) stop("a causal SNP has no variation in the panel")
    X <- sweep(sweep(G, 2, m, "-"), 2, s, "/")
  }
  U <- matrix(NA_real_, n, length(reps)); colnames(U) <- sprintf("rep_%04d", reps)
  zc <- list()
  for (j in seq_along(reps)) {
    r <- reps[j]; set.seed(SEED + 100000L + r)
    u <- stats::rnorm(n) / sqrt(n)
    c_r <- cs[rep_id == r]
    if (nrow(c_r)) {
      ## exact noise z at this replicate's causal SNPs, for step 02's check
      ## (the signal is added in step 02, through local LD only)
      zc[[j]] <- data.table(rep_id = r, snp = c_r$snp,
                            z_exact = as.vector(crossprod(X[, c_r$snp, drop = FALSE], u)))
    }
    U[, j] <- u
  }
  fwrite(cbind(data.table(`#FID` = fam$FID, IID = fam$IID), as.data.table(U)),
         file.path(RUN, "u", sprintf("u_cell%02d.tsv", k)), sep = "\t")
  zc <- if (length(zc)) rbindlist(zc) else data.table(rep_id = integer(0), snp = character(0), z_exact = numeric(0))
  fwrite(zc, file.path(RUN, "u", sprintf("zcheck_cell%02d.tsv", k)), sep = "\t")
  cat(sprintf("[u] cell %d: %d replicates, %d causal SNPs\n", k, length(reps), length(unique(cs$snp))))
}
