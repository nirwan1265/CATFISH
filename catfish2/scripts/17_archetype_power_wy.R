#!/usr/bin/env Rscript
## =============================================================================
## 17_archetype_power_wy.R -- archetype power and omnibus regret, measured under
## the Westfall-Young family-wise decision rule rather than a per-pathway alpha.
##
## WHY THIS IS NOT THE THESIS BLOCK-H ANALYSIS
##   The thesis power analysis (run_block_h in scripts/all_figs.R) simulated ONE
##   pathway in isolation, calibrated it against its own null, and called it a
##   detection when p < 0.05. That is marginal power: the right quantity when a
##   single pathway is tested. This paper reports a family-wise decision over K
##   pathways, so the detection event has to be the family-wise one.
##
##   The archetype layer definitions, the effect-scale mechanism and the
##   permuted gene positions are ported unchanged from all_figs.R so the two
##   analyses differ only in the decision rule. Archetype CODES are relabelled
##   to the manuscript's names (Table 1); the layer parameters are identical:
##       SDA -> SSS   15% at delta 3.00            (strong sparse signal)
##       CME -> MDS   60% at delta 1.10            (moderate dense signal)
##       DPS -> WDE   90% at delta 0.55            (weak diffuse enrichment)
##       HDS -> MSA   10% at 1.80 + 60% at 1.10    (mixed signal architecture)
##       SGP -> SGP    5% at delta 4.00            (single-gene proxy)
##
## DESIGN
##   One shared gene universe with block LD (as in 02/03/04). K pathways of size
##   m are drawn from a common pool, so they overlap. Pathway 1 carries the
##   archetype signal; the other K-1 stay null.
##     calibration  B_CAL draws, everything null, all K pathways
##                  -> per-component and omnibus null banks
##                  -> WY threshold = alpha-quantile of the row-wise MINIMUM
##                     over the K pathways (valid on the statistic scale here
##                     because every pathway has the same size m and LD)
##     evaluation   B_EVAL draws with the signal injected into pathway 1; only
##                  pathway 1 is scored, so evaluation is ~K times cheaper
##   Three detections are recorded for every component and for the omnibus:
##     power_marginal   calibrated p <= alpha            (the thesis quantity)
##     power_wy         statistic   <= WY threshold      (family-wise)
##     power_bonf       analytic  p <  alpha / K      (family-wise competitor)
##   Bonferroni is taken on the ANALYTIC component p, not the calibrated one,
##   exactly as 02_breakage_stress.R::fwer_rates does. A calibrated p from a
##   bank of B draws cannot go below 1/(B+1), which is the same order as
##   alpha/K = 3.3e-4, so a calibrated Bonferroni would be censored at the
##   threshold it is being compared against. The analytic p has no such floor.
##
##   Rscript 17_archetype_power_wy.R [B_CAL] [B_EVAL] [EFFECT_SCALES] [SIZES]
##     EFFECT_SCALES, SIZES: comma-separated, e.g. 1.0,1.5,2.0,2.5
## =============================================================================
source("00_lib.R")
args   <- commandArgs(trailingOnly = TRUE)
B_CAL  <- if (length(args) >= 1) as.integer(args[[1]]) else 2000L
B_EVAL <- if (length(args) >= 2) as.integer(args[[2]]) else 1000L
EFFECTS <- if (length(args) >= 3) as.numeric(strsplit(args[[3]], ",")[[1]]) else
           c(0.5, 1.0, 1.5, 2.0, 2.5, 3.0)
SIZES  <- if (length(args) >= 4) as.integer(strsplit(args[[4]], ",")[[1]]) else
          c(5L, 20L, 50L)
## FALSE drops TFisher from the component set and from the omnibus, so the two
## runs differ only in whether that component is present.
USE_TF <- if (length(args) >= 5) as.logical(args[[5]]) else TRUE
## 6th arg: comma-separated tau grid, overriding 00_lib.R's TAU_GRID_DEFAULT.
## The TFisher authors (Zhang, Tong, Landers & Wu 2020) use TAU1 =
## c(0.01, 0.05, 0.5, 1) and recommend spanning small, medium and large tau;
## the tau=1 corner is dropped here because Fisher is already a component.
if (length(args) >= 6) TAU_GRID_DEFAULT <- as.numeric(strsplit(args[[6]], ",")[[1]])
TAG <- paste0(if (isTRUE(USE_TF)) "" else "_noTF",
              if (length(args) >= 6) sprintf("_tau%s", gsub("[.,]", "", args[[6]])) else "")
set.seed(20260706)

G     <- 3000L     # genes in the universe
BLOCK <- 50L       # LD block size
RHO   <- 0.3       # within-block correlation
K     <- 150L      # pathways in the family
POOL  <- 600L      # pathways are drawn from the first POOL genes
ALPHA <- 0.05
OUT   <- "../results"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

## --- archetypes, ported from all_figs.R::block_a_archetypes ------------------
ARCHETYPES <- list(
  list(code = "SSS", label = "Strong Sparse Signal",
       layers = list(list(prop = 0.15, delta = 3.00))),
  list(code = "MDS", label = "Moderate Dense Signal",
       layers = list(list(prop = 0.60, delta = 1.10))),
  list(code = "WDE", label = "Weak Diffuse Enrichment",
       layers = list(list(prop = 0.90, delta = 0.55))),
  list(code = "MSA", label = "Mixed Signal Architecture",
       layers = list(list(prop = 0.10, delta = 1.80),
                     list(prop = 0.60, delta = 1.10))),
  list(code = "SGP", label = "Single-Gene Proxy",
       layers = list(list(prop = 0.05, delta = 4.00)))
)

## mean shift for one pathway of size m: layer sizes are proportions of m,
## positions permuted per replicate (same as all_figs.R::make_mu_from_layers).
make_mu <- function(m, layers, effect_scale) {
  mu <- rep(0, m)
  avail <- sample.int(m, m)
  start <- 1L
  for (ly in layers) {
    n_ly <- max(1L, as.integer(round(m * ly$prop)))
    n_ly <- min(n_ly, m - start + 1L)
    if (n_ly <= 0L) next
    mu[avail[start:(start + n_ly - 1L)]] <- ly$delta * effect_scale
    start <- start + n_ly
  }
  mu
}

cat(sprintf("archetype power under WY: K=%d pool=%d G=%d B_cal=%d B_eval=%d TFisher=%s\n",
            K, POOL, G, B_CAL, B_EVAL, USE_TF))
cat(sprintf("  sizes   : %s\n", paste(SIZES, collapse = ", ")))
cat(sprintf("  effects : %s\n", paste(EFFECTS, collapse = ", ")))

res <- list()
for (m in SIZES) {
  uni   <- make_universe(G, BLOCK, RHO)
  paths <- make_pathways(K, m, POOL, G)
  sig   <- paths[[1]]                       # the pathway that carries signal
  cat(sprintf("\n  m=%d : calibration bank (%d draws x %d pathways)\n", m, B_CAL, K))

  arr_cal <- components_array(uni, paths, B_CAL, use_tfisher = USE_TF)
  om_cal  <- omnibus_from_components(arr_cal)            # B_CAL x K, log p
  STATS   <- dimnames(arr_cal)[[3]]                      # 5 components + Omnibus
  paths1  <- paths_to_chr(paths, G)[1]                   # the signal pathway

  ## null bank for pathway 1 (marginal calibration) and the WY thresholds
  null_1 <- lapply(STATS, function(s) arr_cal[, 1, s])
  names(null_1) <- STATS
  thr_wy <- vapply(STATS, function(s) wy_threshold(arr_cal[, , s], ALPHA),
                   numeric(1))
  cat("    WY thresholds (log p): ")
  cat(paste(sprintf("%s=%.2f", STATS, thr_wy), collapse = "  "), "\n")

  for (arch in ARCHETYPES) {
    for (eff in EFFECTS) {
      ## evaluation draws: signal in pathway 1 only, so score pathway 1 only
      obs <- matrix(NA_real_, B_EVAL, length(STATS), dimnames = list(NULL, STATS))
      for (b in seq_len(B_EVAL)) {
        z <- draw_universe_z(uni)
        z[sig] <- z[sig] + make_mu(m, arch$layers, eff)
        lp <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
        obs[b, ] <- components_from_draw(lp, z, paths1, G)[1, STATS]
      }

      for (s in STATS) {
        p_marg <- vapply(obs[, s], function(o) mvn_calibrate(o, null_1[[s]]),
                         numeric(1))
        res[[length(res) + 1L]] <- data.frame(
          m = m, archetype = arch$code, archetype_label = arch$label,
          effect_scale = eff, method = s,
          power_marginal = mean(p_marg <= ALPHA, na.rm = TRUE),
          power_wy       = mean(obs[, s] <= thr_wy[[s]], na.rm = TRUE),
          power_bonf     = mean(exp(obs[, s]) <  ALPHA / K, na.rm = TRUE),
          stringsAsFactors = FALSE)
      }
      w <- vapply(STATS, function(s)
             mean(obs[, s] <= thr_wy[[s]], na.rm = TRUE), numeric(1))
      cat(sprintf("    %-4s eff=%.2f  WY power: %s  | best=%s\n",
                  arch$code, eff,
                  paste(sprintf("%s=%.3f", STATS, w), collapse = " "),
                  STATS[which.max(w)]))
    }
  }
}

out <- do.call(rbind, res)
write.csv(out, file.path(OUT, sprintf("17_archetype_power_wy%s.csv", TAG)), row.names = FALSE)

## omnibus regret: how much family-wise power the omnibus gives up against the
## best single component in the same (m, archetype, effect) cell.
cell <- paste(out$m, out$archetype, out$effect_scale)
best <- tapply(out$power_wy[out$method != "Omnibus"],
               cell[out$method != "Omnibus"], max)
om   <- out[out$method == "Omnibus", ]
om$best_component_wy <- as.numeric(best[paste(om$m, om$archetype, om$effect_scale)])
om$regret_wy <- om$best_component_wy - om$power_wy
reg <- do.call(rbind, lapply(split(om, om$archetype), function(d) data.frame(
  archetype = d$archetype[1], archetype_label = d$archetype_label[1],
  mean_regret = mean(d$regret_wy), median_regret = median(d$regret_wy),
  max_regret = max(d$regret_wy),
  within_0.05 = mean(d$regret_wy <= 0.05), within_0.10 = mean(d$regret_wy <= 0.10),
  stringsAsFactors = FALSE)))
write.csv(reg, file.path(OUT, sprintf("17_archetype_regret_wy%s.csv", TAG)), row.names = FALSE)

cat("\nwrote results/17_archetype_power_wy.csv and 17_archetype_regret_wy.csv\n\n")
print(reg, row.names = FALSE)
