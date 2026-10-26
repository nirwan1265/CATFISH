#!/usr/bin/env Rscript
## =============================================================================
## 01_null_calibration.R -- type I error and genomic-control lambda under the
## global null, across LD strength and pathway size, for every calibration.
##
##   analytic    : component p-values combined directly (assumes independence)
##   mvn         : per-pathway MVN calibration (CATFISH's current default)
##   bonferroni  : analytic, threshold alpha/NP
##   wy          : genome-wide max-statistic (this paper)
##
##   Rscript 01_null_calibration.R [NREP] [B] [NP]
## =============================================================================
source("00_lib.R")
args  <- commandArgs(trailingOnly = TRUE)
NREP  <- if (length(args) >= 1) as.integer(args[[1]]) else 500L   # null datasets
B     <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L  # calibration draws
NP    <- if (length(args) >= 3) as.integer(args[[3]]) else 200L   # pathways per dataset
set.seed(20260706)

RHO  <- c(independent = 0.0, moderate = 0.3, strong = 0.6)
SIZE <- c(5L, 20L, 50L)
ALPHA <- c(0.05, 0.01)
OUT <- "../results"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

cat(sprintf("null calibration: %d datasets x %d pathways, B=%d\n", NREP, NP, B))
res <- list()

for (rn in names(RHO)) {
  for (m in SIZE) {
    U <- safe_chol(R_exch(m, RHO[[rn]]))
    ## shared per-size null bank for MVN calibration
    bank <- vapply(seq_len(B), function(i) omnibus_log(draw_lp(U)), numeric(1))

    p_analytic <- numeric(NREP * NP)
    p_mvn      <- numeric(NREP * NP)
    minp_wy    <- numeric(NREP)
    k <- 0L
    for (r in seq_len(NREP)) {
      st <- vapply(seq_len(NP), function(j) omnibus_log(draw_lp(U)), numeric(1))
      idx <- (k + 1L):(k + NP); k <- k + NP
      p_analytic[idx] <- exp(st)
      p_mvn[idx] <- vapply(st, function(s) mvn_calibrate(s, bank), numeric(1))
      minp_wy[r] <- min(st)               # max-statistic across the family
    }
    ## WY: family-wise threshold from the distribution of the minimum
    wy_thresh <- stats::quantile(minp_wy, ALPHA, names = FALSE)

    for (ai in seq_along(ALPHA)) {
      a <- ALPHA[ai]
      res[[length(res) + 1L]] <- data.frame(
        rho = RHO[[rn]], ld = rn, m = m, alpha = a,
        t1_analytic   = type1(p_analytic, a),
        t1_mvn        = type1(p_mvn, a),
        t1_bonferroni = type1(p_analytic, a / NP) * NP,  # family-wise, rescaled
        t1_wy         = mean(minp_wy <= wy_thresh[ai]),
        lambda_analytic = lambda_gc(p_analytic),
        lambda_mvn      = lambda_gc(p_mvn),
        stringsAsFactors = FALSE)
    }
    cat(sprintf("  LD=%-11s m=%-3d  analytic t1(.05)=%.4f  mvn=%.4f  lambda_an=%.2f\n",
                rn, m, type1(p_analytic, 0.05), type1(p_mvn, 0.05), lambda_gc(p_analytic)))
  }
}
out <- do.call(rbind, res)
write.csv(out, file.path(OUT, "01_null_calibration.csv"), row.names = FALSE)
cat("\nwrote results/01_null_calibration.csv\n")
print(out[out$alpha == 0.05, c("ld","m","t1_analytic","t1_mvn","t1_wy","lambda_analytic","lambda_mvn")])
