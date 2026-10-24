suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(readr)
  library(scales)
})

set.seed(20260803)

out_dir_fig <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/supp_figs"
out_dir_tab <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_sup/tables"
dir.create(out_dir_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tab, recursive = TRUE, showWarnings = FALSE)

n_reps <- 400L
n_genes <- 500L
n_pathways <- 50L
pathway_size_range <- 1:15
n_top_pathways <- 10L
n_true_pathways <- 5L
n_causal_in_path <- 30L
n_causal_background <- 5L
top_k <- 20L

w_path_grid <- c(0, 0.25, 0.5, 0.75, 1.0)
w_hub_grid  <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5)

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
      best_path_p[g] <- min(p_path[containsing <- containing])
      top_path_count[g] <- length(containing)
    }
  }

  s_gwas <- rank_score(p_gwas)
  s_magma <- rank_score(p_magma)
  s_path <- rank_score(best_path_p, supported = !is.na(best_path_p))
  b_hub <- log10(top_path_count + 1)

  tibble(
    gene = genes,
    causal = causal,
    s_gwas = s_gwas,
    s_magma = s_magma,
    s_path = s_path,
    b_hub = b_hub
  )
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
  tibble(
    w_path = w_path,
    w_hub = w_hub,
    baseline_overlap = baseline_overlap,
    top20_causal_fraction = top_causal_fraction,
    rank_correlation = rank_corr
  )
}

message("Running candidate-gene weight sensitivity simulation...")

all_results <- vector("list", n_reps)

for (rep_id in seq_len(n_reps)) {
  dat <- simulate_one_dataset()
  baseline <- score_and_rank(dat, 0.5, 0.3)
  baseline_top <- baseline$gene[seq_len(top_k)]

  grid_res <- bind_rows(lapply(w_path_grid, function(wp) {
    bind_rows(lapply(w_hub_grid, function(wh) {
      evaluate_weights(dat, wp, wh, baseline_top)
    }))
  })) %>%
    mutate(rep = rep_id)

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
    .groups = "drop"
  )

readr::write_csv(results_long,
                 file.path(out_dir_tab, "gene_priority_weight_sensitivity_replicates.csv"))
readr::write_csv(summary_res,
                 file.path(out_dir_tab, "gene_priority_weight_sensitivity_summary.csv"))

make_heat <- function(df, fill_col, label_col, title, low = "#f7fbff", high = "#08306b") {
  ggplot(df, aes(x = factor(w_path), y = factor(w_hub), fill = .data[[fill_col]])) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = sprintf("%.2f", .data[[label_col]])), size = 4.5, fontface = "bold") +
    scale_fill_gradient(low = low, high = high, limits = c(min(df[[fill_col]]), max(df[[fill_col]]))) +
    labs(
      title = title,
      x = "Pathway weight",
      y = "Hub weight",
      fill = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    )
}

p_overlap <- make_heat(
  summary_res,
  fill_col = "mean_top20_overlap",
  label_col = "mean_top20_overlap",
  title = "Top-20 overlap with baseline weights (0.5, 0.3)",
  low = "#fff5f0",
  high = "#cb181d"
)

p_causal <- make_heat(
  summary_res,
  fill_col = "mean_top20_causal_fraction",
  label_col = "mean_top20_causal_fraction",
  title = "Fraction of causal genes recovered in top 20",
  low = "#f7fcf5",
  high = "#238b45"
)

p_corr <- make_heat(
  summary_res,
  fill_col = "mean_rank_correlation",
  label_col = "mean_rank_correlation",
  title = "Spearman correlation with baseline full ranking",
  low = "#f7fbff",
  high = "#08519c"
)

fig <- (p_overlap | p_causal | p_corr) +
  plot_annotation(
    title = "Candidate-gene prioritization weight sensitivity (truth-planted simulation)",
    subtitle = paste0(
      n_reps, " replicates; 500 genes, 50 pathways, 5 enriched pathways, top-", top_k,
      " baseline set defined at pathway weight 0.5 and hub weight 0.3"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    )
  )

ggsave(
  filename = file.path(out_dir_fig, "gene_priority_weight_sensitivity_simulation.png"),
  plot = fig,
  width = 17,
  height = 6.5,
  dpi = 350,
  bg = "white"
)

message("Done.")
print(summary_res %>% arrange(desc(mean_top20_overlap)))
