#!/usr/bin/env Rscript
## =============================================================================
## 22_weight_sensitivity.R -- sensitivity of the candidate-gene score to its
## pathway and hub weights (Supplementary Figure, Section 2.5).
##
##   score = S_GWAS + S_MAGMA + w_path * S_PATH + w_hub * B_hub
##   grid:  w_path in {0, 0.25, 0.5, 0.75, 1},  w_hub in {0, 0.1, ..., 0.5}
##
## Part A  truth-planted simulation. The simulation code, seed and settings are
##         those of Final/final_scripts/gene_priority_weight_sensitivity_sim.R
##         (thesis Appendix A), so it reproduces the thesis numbers exactly.
##         400 replicates of 500 genes and 50 pathways with planted causal genes.
## Part B  the real lipid analyses: the same grid applied to the per-layer
##         scores already computed by 06_candidate_genes.R, summarised by the
##         ranks of the pre-specified established genes. This is a sensitivity
##         report, not a tuning step: the weights were fixed beforehand.
##
## Run from catfish2/scripts after 06 (SPECIES=human, TRAIT=ldl and hdl):
##   Rscript 22_weight_sensitivity.R
##   N_REPS=50 Rscript 22_weight_sensitivity.R      # quick check
## Writes ../results/22_weight_sensitivity_sim_{replicates,summary}.csv
##        ../results/22_weight_sensitivity_lipids.csv
##        ../figures/Fig16_weight_sensitivity.png
## =============================================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
})

RES <- "../results"; FIG <- "../figures"
dir.create(RES, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)

w_path_grid <- c(0, 0.25, 0.5, 0.75, 1.0)
w_hub_grid  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)
W_PATH0 <- 0.5; W_HUB0 <- 0.3     # the defaults used throughout

## ======================================================= Part A: simulation ==
## --- verbatim from gene_priority_weight_sensitivity_sim.R (thesis) ---------
set.seed(20260803)

n_reps <- as.integer(Sys.getenv("N_REPS", "400"))
n_genes <- 500L
n_pathways <- 50L
pathway_size_range <- 1:15
n_top_pathways <- 10L
n_true_pathways <- 5L
n_causal_in_path <- 30L
n_causal_background <- 5L
top_k <- 20L

rank_score <- function(p, supported = rep(TRUE, length(p))) {
  out <- numeric(length(p))
  idx <- which(supported & is.finite(p) & !is.na(p))
  if (length(idx) == 0) return(out)
  r <- rank(p[idx], ties.method = "average")
  out[idx] <- -log10(r / length(idx))
  out
}

simulate_one_dataset <- function() {
  genes <- paste0("G", seq_len(n_genes))

  pathway_sizes <- sample(pathway_size_range, n_pathways, replace = TRUE)
  pathways <- lapply(pathway_sizes, function(sz) sample.int(n_genes, sz, replace = FALSE))
  names(pathways) <- paste0("P", seq_len(n_pathways))

  eligible_true <- which(pathway_sizes >= 6)
  if (length(eligible_true) < n_true_pathways) eligible_true <- seq_len(n_pathways)
  true_path_idx <- sample(eligible_true, n_true_pathways, replace = FALSE)

  true_union <- sort(unique(unlist(pathways[true_path_idx])))
  membership_true <- integer(n_genes)
  for (j in true_path_idx) {
    membership_true[pathways[[j]]] <- membership_true[pathways[[j]]] + 1L
  }

  prob_true <- membership_true[true_union] + 0.25
  prob_true <- prob_true / sum(prob_true)
  causal_in_path <- sample(true_union,
                           size = min(n_causal_in_path, length(true_union)),
                           replace = FALSE,
                           prob = prob_true)

  background_pool <- setdiff(seq_len(n_genes), true_union)
  causal_background <- if (length(background_pool) > 0) {
    sample(background_pool,
           size = min(n_causal_background, length(background_pool)),
           replace = FALSE)
  } else integer(0)

  causal_idx <- sort(unique(c(causal_in_path, causal_background)))
  causal <- integer(n_genes)
  causal[causal_idx] <- 1L

  shared_noise <- rnorm(n_genes)
  z_gwas <- 1.9 * causal + 0.20 * membership_true + shared_noise + rnorm(n_genes, sd = 0.8)
  z_magma <- 1.4 * causal + 0.35 * membership_true + 0.6 * shared_noise + rnorm(n_genes, sd = 0.7)

  p_gwas <- pnorm(z_gwas, lower.tail = FALSE)
  p_magma <- pnorm(z_magma, lower.tail = FALSE)

  path_z <- numeric(n_pathways)
  for (j in seq_len(n_pathways)) {
    members <- pathways[[j]]
    c_count <- sum(causal[members])
    c_frac <- mean(causal[members])
    in_true_set <- as.integer(j %in% true_path_idx)
    path_z[j] <- 1.2 * in_true_set + 0.35 * c_count + 1.1 * c_frac + rnorm(1, sd = 0.9)
  }
  p_path <- pnorm(path_z, lower.tail = FALSE)

  top_path_idx <- order(p_path)[seq_len(min(n_top_pathways, n_pathways))]

  best_path_p <- rep(NA_real_, n_genes)
  top_path_count <- integer(n_genes)
  for (g in seq_len(n_genes)) {
    containing <- top_path_idx[vapply(pathways[top_path_idx], function(m) g %in% m, logical(1))]
    if (length(containing) > 0) {
      best_path_p[g] <- min(p_path[containing])
      top_path_count[g] <- length(containing)
    }
  }

  s_gwas <- rank_score(p_gwas)
  s_magma <- rank_score(p_magma)
  s_path <- rank_score(best_path_p, supported = !is.na(best_path_p))
  b_hub <- log10(top_path_count + 1)

  tibble(gene = genes, causal = causal, s_gwas = s_gwas, s_magma = s_magma,
         s_path = s_path, b_hub = b_hub)
}

score_and_rank <- function(dat, w_path, w_hub) {
  scored <- dat %>%
    mutate(score = s_gwas + s_magma + w_path * s_path + w_hub * b_hub) %>%
    arrange(desc(score), desc(s_magma), desc(s_gwas))
  scored$rank <- seq_len(nrow(scored))
  scored
}

evaluate_weights <- function(dat, w_path, w_hub, baseline_top) {
  scored <- score_and_rank(dat, w_path, w_hub)
  top_genes <- scored$gene[seq_len(top_k)]
  baseline_overlap <- length(intersect(top_genes, baseline_top)) / top_k
  top_causal_fraction <- mean(scored$causal[seq_len(top_k)])
  rank_corr <- suppressWarnings(cor(scored$rank,
                                    score_and_rank(dat, 0.5, 0.3)$rank[match(scored$gene,
                                                                             score_and_rank(dat, 0.5, 0.3)$gene)],
                                    method = "spearman"))
  tibble(w_path = w_path, w_hub = w_hub, baseline_overlap = baseline_overlap,
         top20_causal_fraction = top_causal_fraction, rank_correlation = rank_corr)
}

message(sprintf("Part A: %d simulation replicates ...", n_reps))
all_results <- vector("list", n_reps)
for (rep_id in seq_len(n_reps)) {
  dat <- simulate_one_dataset()
  baseline <- score_and_rank(dat, 0.5, 0.3)
  baseline_top <- baseline$gene[seq_len(top_k)]
  grid_res <- bind_rows(lapply(w_path_grid, function(wp) {
    bind_rows(lapply(w_hub_grid, function(wh) evaluate_weights(dat, wp, wh, baseline_top)))
  })) %>% mutate(rep = rep_id)
  all_results[[rep_id]] <- grid_res
}
results_long <- bind_rows(all_results)
summary_res <- results_long %>%
  group_by(w_path, w_hub) %>%
  summarise(
    mean_top20_overlap = mean(baseline_overlap, na.rm = TRUE),
    se_top20_overlap = sd(baseline_overlap, na.rm = TRUE) / sqrt(n()),
    mean_top20_causal_fraction = mean(top20_causal_fraction, na.rm = TRUE),
    se_top20_causal_fraction = sd(top20_causal_fraction, na.rm = TRUE) / sqrt(n()),
    mean_rank_correlation = mean(rank_correlation, na.rm = TRUE),
    se_rank_correlation = sd(rank_correlation, na.rm = TRUE) / sqrt(n()),
    .groups = "drop")
## --- end verbatim -----------------------------------------------------------
utils::write.csv(results_long, file.path(RES, "22_weight_sensitivity_sim_replicates.csv"), row.names = FALSE)
utils::write.csv(summary_res,  file.path(RES, "22_weight_sensitivity_sim_summary.csv"),    row.names = FALSE)

## ============================================ Part B: the lipid analyses ====
CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))
lip <- list()
for (t in names(CONTROLS)) {
  f <- file.path(RES, sprintf("06_candidate_genes_%s.csv", t))
  if (!file.exists(f)) { message("Part B: ", f, " not found; skipping ", t); next }
  g <- utils::read.csv(f, stringsAsFactors = FALSE)
  ci <- which(g$symbol %in% CONTROLS[[t]])
  for (wp in w_path_grid) for (wh in w_hub_grid) {
    sc <- g$S_GWAS + g$S_MAGMA + wp * g$S_PATH + wh * g$B_hub
    rk <- rank(-sc, ties.method = "first")
    r  <- rk[ci]
    lip[[length(lip) + 1L]] <- data.frame(trait = t, w_path = wp, w_hub = wh,
      n_controls = length(ci), median_rank = stats::median(r),
      top10 = sum(r <= 10), top50 = sum(r <= 50))
  }
}
lip <- do.call(rbind, lip)
if (!is.null(lip))
  utils::write.csv(lip, file.path(RES, "22_weight_sensitivity_lipids.csv"), row.names = FALSE)

## ================================================================ figure ====
heat <- function(z, main, fmt, lower_is_better = FALSE) {
  ## z: matrix [w_path x w_hub]
  pal <- grDevices::colorRampPalette(c("#f7fbff", "#6baed6", "#08306b"))(64)
  v <- if (lower_is_better) -z else z
  rng <- range(v); if (diff(rng) == 0) rng <- rng + c(-1, 1)
  col <- pal[1 + round(63 * (v - rng[1]) / diff(rng))]
  nr <- nrow(z); nc <- ncol(z)
  plot(NA, xlim = c(0.5, nr + 0.5), ylim = c(0.5, nc + 0.5), xaxt = "n", yaxt = "n",
       xlab = "pathway weight", ylab = "hub weight", main = main, xaxs = "i", yaxs = "i")
  for (i in seq_len(nr)) for (j in seq_len(nc)) {
    rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col = col[(j - 1) * nr + i], border = "white")
    vv <- v[i, j]; dark <- (vv - rng[1]) / diff(rng) > 0.55
    text(i, j, sprintf(fmt, z[i, j]), cex = 0.62, col = if (dark) "white" else "black")
  }
  i0 <- match(W_PATH0, w_path_grid); j0 <- match(W_HUB0, w_hub_grid)
  rect(i0 - 0.5, j0 - 0.5, i0 + 0.5, j0 + 0.5, border = "#C44E52", lwd = 2.5)
  axis(1, at = seq_len(nr), labels = w_path_grid, cex.axis = 0.75)
  axis(2, at = seq_len(nc), labels = w_hub_grid, las = 1, cex.axis = 0.75)
}
as_mat <- function(d, col) {
  m <- matrix(NA_real_, length(w_path_grid), length(w_hub_grid))
  m[cbind(match(d$w_path, w_path_grid), match(d$w_hub, w_hub_grid))] <- d[[col]]
  m
}
s <- as.data.frame(summary_res)
npan <- 2L + if (is.null(lip)) 0L else length(unique(lip$trait))
png(file.path(FIG, "Fig16_weight_sensitivity.png"), width = 600 * npan, height = 620, res = 165)
op <- par(mfrow = c(1, npan), mar = c(4, 4, 3.2, 0.8), mgp = c(2.3, 0.6, 0))
heat(as_mat(s, "mean_top20_causal_fraction"),
     "A  simulation: causal genes in top 20", "%.2f")
heat(as_mat(s, "mean_rank_correlation"),
     "B  simulation: rank correlation\nwith the default weights", "%.3f")
lab <- c(ldl = "LDL-C", hdl = "HDL-C"); k <- 2L
if (!is.null(lip)) for (t in unique(lip$trait)) {
  k <- k + 1L
  heat(as_mat(lip[lip$trait == t, ], "median_rank"),
       sprintf("%s  %s: median rank of\n%d established genes", LETTERS[k], lab[[t]],
               lip$n_controls[lip$trait == t][1]), "%g", lower_is_better = TRUE)
}
par(op); dev.off()
cat("wrote ../figures/Fig16_weight_sensitivity.png\n")

## ---------------------------------------------------------- console summary --
d0 <- s[s$w_path == W_PATH0 & s$w_hub == W_HUB0, ]
cat(sprintf("\nsimulation (%d replicates): causal fraction in top 20 %.3f at defaults; range %.3f-%.3f; min rank correlation with defaults %.3f; min top-20 overlap %.3f\n",
            n_reps, d0$mean_top20_causal_fraction, min(s$mean_top20_causal_fraction),
            max(s$mean_top20_causal_fraction), min(s$mean_rank_correlation),
            min(s$mean_top20_overlap)))
if (!is.null(lip)) print(lip[order(lip$trait, lip$w_path, lip$w_hub), ], row.names = FALSE)
