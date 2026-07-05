#!/usr/bin/env Rscript
# =============================================================================
# 08_aggregate_and_plot.R  --  Pool null pathway p-values across all runs
# =============================================================================
suppressMessages({ library(data.table); library(ggplot2) })

args <- commandArgs(trailingOnly = TRUE)
in_dir  <- if (length(args) >= 1) args[1] else Sys.getenv("CATFISH_OUT")
out_dir <- if (length(args) >= 2) args[2] else Sys.getenv("FINAL_DIR")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

files <- list.files(in_dir, pattern = "^pathway_pvals_perm_.*\\.csv$", full.names = TRUE)
if (length(files) == 0) stop(paste("no per-perm CSVs in", in_dir))
cat(sprintf("[agg] pooling %d permutation files\n", length(files)))
dat <- rbindlist(lapply(files, fread), fill = TRUE)

pcols <- intersect(c("omni_p_final","omni_p_mvn","omni_p_analytic",
                     "acat_p","fisher_p","minp_p_analytic",
                     "stouffer_p_analytic","tfisher_p_analytic",
                     "acat_p_mvn_cal","fisher_p_mvn_cal","tfisher_p_mvn_cal",
                     "minp_p_mvn_cal","stouffer_p_mvn_cal",
                     "omni_p_mvn_compcal"), names(dat))

type1 <- rbindlist(lapply(pcols, function(pc) {
  p <- dat[[pc]]; p <- p[is.finite(p) & p >= 0 & p <= 1]
  n <- length(p)
  if (n < 1L) return(NULL)
  one <- function(a) {
    k <- sum(p <= a); rate <- k / n
    ci <- binom.test(k, n, a)$conf.int
    data.table(method = pc, alpha = a, n = n, observed = rate,
               ci_lo = ci[1], ci_hi = ci[2],
               calibrated = (a >= ci[1] & a <= ci[2]))
  }
  rbind(one(0.05), one(0.01))
}))
fwrite(type1, file.path(out_dir, "type1_error_table.csv"))

lambda_of <- function(p){ p<-p[is.finite(p)&p>0&p<=1]
  if (!length(p)) return(NA_real_)
  median(qchisq(1-p,1))/qchisq(0.5,1) }
lam <- data.table(method = pcols, lambda = sapply(pcols, function(pc) lambda_of(dat[[pc]])))
lam <- lam[is.finite(lambda) & !is.na(lambda)]
fwrite(lam, file.path(out_dir, "lambda_table.csv"))

qq_dt <- function(p){ p<-sort(p[is.finite(p)&p>0&p<=1]); n<-length(p)
  data.table(expected=-log10(ppoints(n)), observed=-log10(p)) }
head_p <- if ("omni_p_final" %in% pcols) "omni_p_final" else pcols[1]
qd <- qq_dt(dat[[head_p]])
ci_band <- {
  n <- nrow(qd); k <- seq_len(n)
  data.table(expected = -log10(ppoints(n)),
             lo = -log10(qbeta(0.975, k, n - k + 1)),
             hi = -log10(qbeta(0.025, k, n - k + 1)))
}
p_qq <- ggplot() +
  geom_ribbon(data = ci_band, aes(expected, ymin = lo, ymax = hi),
              fill = "grey80", alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, colour = "red", linetype = "dashed") +
  geom_point(data = qd, aes(expected, observed), size = 0.7, alpha = 0.5) +
  labs(x = expression(Expected~-log[10](p)),
       y = expression(Observed~-log[10](p)),
       title = "CATFISH omnibus under spike simulation",
       subtitle = sprintf("%d simulations, %d pooled pathway p-values; lambda = %.3f",
                          length(files), nrow(qd), lambda_of(dat[[head_p]]))) +
  theme_minimal(base_size = 12) + theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "QQ_omnibus_null.png"), p_qq, width = 6, height = 6, dpi = 200)
ggsave(file.path(out_dir, "QQ_omnibus_null.pdf"), p_qq, width = 6, height = 6)

p_hist <- ggplot(data.table(p = dat[[head_p]][is.finite(dat[[head_p]])]),
                 aes(p)) +
  geom_histogram(breaks = seq(0,1,0.05), fill = "steelblue", colour = "white") +
  geom_hline(yintercept = sum(is.finite(dat[[head_p]]))/20,
             colour = "red", linetype = "dashed") +
  labs(x = "Omnibus p-value", y = "Count",
       title = "P-value histogram") +
  theme_minimal(base_size = 12)
ggsave(file.path(out_dir, "hist_omnibus_null.png"), p_hist, width = 6, height = 4, dpi = 200)

cat(sprintf("[agg] figures + tables written to %s\n", out_dir))
