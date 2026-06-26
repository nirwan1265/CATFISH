suppressPackageStartupMessages({
  source("/Users/nirwantandukar/Documents/Github/MAGCAT/scripts/all_figs.R")
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(stringr)
  library(scales)
})

# ------------------------------------------------------------------------------
# Final Figure 7 script
# Dry tons per acre: CATFISH method comparison + integrated candidate evidence
# Layout:
#   A. UpSet-style overlap among significant CATFISH method hits
#   B. Omnibus QQ plot with final pathway significance categories
#   C. Method concordance matrix (Spearman upper, Jaccard lower)
#   D. Integrated GWAS-MAGMA-pathway candidate-gene prioritization
# ------------------------------------------------------------------------------

final_fig_dir <- "/Users/nirwantandukar/Documents/Github/MAGCAT/Final/final_main_fig"
dir.create(final_fig_dir, recursive = TRUE, showWarnings = FALSE)

base_dir <- "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/Dry_tons_per_acre"
catfish_file <- file.path(
  base_dir,
  "CATFISH_permutation_B10000_mvn_GPD_paper_tau_false",
  "Dry_tons_per_acre_CATFISH_ACAT_mvn_B10000_GPD_paper_tau_false.csv"
)
gene_score_file <- file.path(
  base_dir,
  "candidate_gene_scoring_B10000_GPD_paper_tau_false",
  "candidate_genes_all_Dry_tons_per_acre_B10000_GPD_paper_tau_false.csv"
)

if (!file.exists(catfish_file)) stop("CATFISH file not found: ", catfish_file, call. = FALSE)
if (!file.exists(gene_score_file)) stop("Candidate gene file not found: ", gene_score_file, call. = FALSE)

panel_theme <- plot_theme +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 14, face = "bold"),
    axis.text.x = element_text(size = 11, color = "black"),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.line = element_line(color = "black", linewidth = 0.35),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    plot.margin = margin(6, 8, 6, 8)
  )

pick_component_col <- function(df, candidates) {
  for (col in candidates) {
    if (!col %in% names(df)) next
    vals <- suppressWarnings(as.numeric(df[[col]]))
    if (any(is.finite(vals) & !is.na(vals))) return(col)
  }
  stop(
    "Could not find a populated component column among: ",
    paste(candidates, collapse = ", "),
    call. = FALSE
  )
}

jaccard_index <- function(a, b) {
  if (!length(a) && !length(b)) return(1)
  u <- union(a, b)
  if (!length(u)) return(0)
  length(intersect(a, b)) / length(u)
}

catfish_df <- fread(catfish_file, data.table = FALSE)
gene_df <- read.csv(gene_score_file, stringsAsFactors = FALSE, check.names = FALSE)

component_cols <- c(
  ACAT = pick_component_col(catfish_df, c("acat_p", "acat_p_mvn_cal")),
  Fisher = pick_component_col(catfish_df, c("fisher_p", "fisher_p_mvn_cal")),
  TFisher = pick_component_col(catfish_df, c("tfisher_p_analytic", "tfisher_p_mvn_cal")),
  minP = pick_component_col(catfish_df, c("minp_p_analytic", "minp_p_mvn_cal")),
  Stouffer = pick_component_col(catfish_df, c("stouffer_p_analytic", "stouffer_p_mvn_cal")),
  Omnibus = pick_component_col(catfish_df, c("omni_p_final", "omni_p_mvn", "omni_p_analytic"))
)

method_order <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")
upset_method_order <- c("Stouffer", "minP", "TFisher", "Fisher", "ACAT", "Omnibus")

n_pathways <- nrow(catfish_df)
bonf_thresh <- 0.05 / n_pathways

sig_sets <- lapply(method_order, function(method) {
  col <- component_cols[[method]]
  vals <- suppressWarnings(as.numeric(catfish_df[[col]]))
  catfish_df$pathway_id[is.finite(vals) & !is.na(vals) & vals < bonf_thresh]
})
names(sig_sets) <- method_order

# ------------------------------------------------------------------------------
# Panel A: UpSet-style overlap plot
# ------------------------------------------------------------------------------
membership_df <- data.frame(pathway_id = catfish_df$pathway_id, stringsAsFactors = FALSE)
for (method in method_order) {
  membership_df[[method]] <- membership_df$pathway_id %in% sig_sets[[method]]
}
membership_df <- membership_df %>%
  filter(if_any(all_of(method_order), identity))

combo_df <- membership_df %>%
  mutate(combo_key = apply(dplyr::select(., all_of(method_order)), 1, function(x) paste(as.integer(x), collapse = ""))) %>%
  count(combo_key, sort = TRUE, name = "intersection_size") %>%
  mutate(
    combo_id = row_number(),
    active_methods = lapply(combo_key, function(key) method_order[as.logical(as.integer(strsplit(key, "")[[1]]))])
  ) %>%
  slice_head(n = 8)

set_sizes_df <- data.frame(
  method = factor(upset_method_order, levels = upset_method_order),
  set_size = sapply(upset_method_order, function(m) sum(membership_df[[m]], na.rm = TRUE)),
  fill_group = ifelse(upset_method_order == "Omnibus", "Omnibus", "Component"),
  stringsAsFactors = FALSE
)

combo_matrix_df <- expand.grid(
  combo_id = combo_df$combo_id,
  method = upset_method_order,
  stringsAsFactors = FALSE
) %>%
  left_join(combo_df %>% dplyr::select(combo_id, active_methods), by = "combo_id") %>%
  mutate(
    method = factor(method, levels = upset_method_order),
    y_pos = match(method, upset_method_order),
    active = mapply(function(m, a) m %in% a, method, active_methods),
    combo_id = factor(combo_id, levels = combo_df$combo_id)
  )

combo_segments_df <- combo_matrix_df %>%
  filter(active) %>%
  group_by(combo_id) %>%
  summarise(
    y_min = min(as.integer(method)),
    y_max = max(as.integer(method)),
    .groups = "drop"
  )

row_bg_df <- data.frame(
  ymin = seq_along(upset_method_order) - 0.5,
  ymax = seq_along(upset_method_order) + 0.5,
  fill = rep(c("odd", "even"), length.out = length(upset_method_order))
)

panel_a_top <- ggplot(combo_df, aes(x = factor(combo_id, levels = combo_id), y = intersection_size)) +
  geom_col(fill = "#3F3F3F", width = 0.62) +
  geom_text(aes(label = intersection_size), vjust = -0.4, size = 3.4, color = "#3F3F3F") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.10))) +
  labs(title = NULL, x = NULL, y = "Intersection Size") +
  panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  )

panel_a_sets <- ggplot(set_sizes_df, aes(x = set_size, y = method, fill = fill_group)) +
  geom_col(width = 0.62) +
  scale_fill_manual(values = c(Component = "#3E6FA5", Omnibus = "#A31818")) +
  scale_y_discrete(position = "right") +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Set Size", y = NULL) +
  panel_theme +
  theme(
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 11, color = "black"),
    axis.text.x = element_text(size = 10),
    legend.position = "none",
    plot.title = element_blank()
  )

panel_a_matrix <- ggplot() +
  geom_rect(
    data = row_bg_df,
    aes(xmin = 0.5, xmax = length(combo_df$combo_id) + 0.5, ymin = ymin, ymax = ymax, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.18
  ) +
  scale_fill_manual(values = c(odd = "#F3F3F3", even = "#FFFFFF")) +
  geom_segment(
    data = combo_segments_df,
    aes(x = combo_id, xend = combo_id, y = y_min, yend = y_max),
    linewidth = 0.7,
    color = "#333333"
  ) +
  geom_point(
    data = combo_matrix_df,
    aes(x = combo_id, y = y_pos),
    color = "#D9D9D9",
    size = 2.8
  ) +
  geom_point(
    data = combo_matrix_df %>% filter(active),
    aes(x = combo_id, y = y_pos),
    color = "#333333",
    size = 3.2
  ) +
  scale_y_continuous(
    breaks = seq_along(upset_method_order),
    labels = rep("", length(upset_method_order)),
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  labs(x = NULL, y = NULL) +
  panel_theme +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none",
    plot.title = element_blank()
  )

panel_a_inner <- (plot_spacer() | panel_a_top) / (panel_a_sets | panel_a_matrix) +
  plot_layout(widths = c(1.1, 2.5), heights = c(1.5, 1), guides = "collect")

# ------------------------------------------------------------------------------
# Panel B: Ranked pathway signal profiles across methods
# ------------------------------------------------------------------------------
panel_b_cols <- c(
  ACAT = "acat_p",
  Fisher = "fisher_p",
  TFisher = "tfisher_p_analytic",
  minP = "minp_p_analytic",
  Stouffer = "stouffer_p_analytic",
  Omnibus = "omni_p_final"
)

missing_panel_b <- setdiff(unname(panel_b_cols), names(catfish_df))
if (length(missing_panel_b)) {
  stop(
    "Panel B columns not found in CATFISH results: ",
    paste(missing_panel_b, collapse = ", "),
    call. = FALSE
  )
}

panel_b_df <- bind_rows(lapply(names(panel_b_cols), function(method) {
  col <- panel_b_cols[[method]]
  pvals <- suppressWarnings(as.numeric(catfish_df[[col]]))
  keep <- is.finite(pvals) & !is.na(pvals) & pvals > 0
  df <- data.frame(
    pathway_id = catfish_df$pathway_id[keep],
    method = method,
    p_value = pvals[keep],
    stringsAsFactors = FALSE
  )
  df <- df[order(df$p_value), , drop = FALSE]
  df$rank <- seq_len(nrow(df))
  df$logP <- -log10(df$p_value)
  df$sig_class <- ifelse(df$p_value < bonf_thresh, "Bonferroni < 0.05", "Not significant")
  df
}))

panel_b_df$method <- factor(panel_b_df$method, levels = names(panel_b_cols))

panel_b_colors <- c(
  ACAT = "#D55E00",
  Fisher = "#CC79A7",
  TFisher = "#7B61FF",
  minP = "#E6AB02",
  Stouffer = "#1B9E77",
  Omnibus = "#B22222",
  "Not significant" = "#B8B8B8"
)

panel_b <- ggplot(panel_b_df, aes(x = rank, y = logP, color = method)) +
  geom_line(linewidth = 1, alpha = 0.9) +
  scale_color_manual(values = panel_b_colors[names(panel_b_cols)], drop = FALSE) +
  labs(
    title = NULL,
    x = "Pathway Rank",
    y = expression(-log[10](italic(p))),
    caption = paste0(
      "n = ", n_pathways,
      " pathways per method"
    )
  ) +
  scale_x_continuous(
    limits = c(1, n_pathways),
    breaks = pretty_breaks(n = 5)
  ) +
  scale_y_continuous(
    limits = c(0, max(panel_b_df$logP, na.rm = TRUE) * 1.05),
    breaks = pretty_breaks(n = 5)
  ) +
  panel_theme +
  theme(
    legend.position = "top",
    plot.caption = element_text(size = 9, color = "#444444")
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 3, alpha = 1)
    )
  )

# ------------------------------------------------------------------------------
# Panel C: Method concordance matrix
# ------------------------------------------------------------------------------
method_p_df <- catfish_df %>%
  transmute(
    pathway_id = pathway_id,
    ACAT = suppressWarnings(as.numeric(.data[[component_cols[["ACAT"]]]])),
    Fisher = suppressWarnings(as.numeric(.data[[component_cols[["Fisher"]]]])),
    TFisher = suppressWarnings(as.numeric(.data[[component_cols[["TFisher"]]]])),
    minP = suppressWarnings(as.numeric(.data[[component_cols[["minP"]]]])),
    Stouffer = suppressWarnings(as.numeric(.data[[component_cols[["Stouffer"]]]])),
    Omnibus = suppressWarnings(as.numeric(.data[[component_cols[["Omnibus"]]]]))
  )

spearman_mat <- matrix(NA_real_, nrow = length(method_order), ncol = length(method_order),
                       dimnames = list(method_order, method_order))
jaccard_mat <- spearman_mat

for (i in seq_along(method_order)) {
  for (j in seq_along(method_order)) {
    mi <- method_order[[i]]
    mj <- method_order[[j]]
    xi <- -log10(pmax(method_p_df[[mi]], 1e-50))
    xj <- -log10(pmax(method_p_df[[mj]], 1e-50))
    spearman_mat[i, j] <- suppressWarnings(cor(xi, xj, method = "spearman", use = "pairwise.complete.obs"))
    jaccard_mat[i, j] <- jaccard_index(sig_sets[[mi]], sig_sets[[mj]])
  }
}

matrix_long <- expand.grid(
  MethodX = method_order,
  MethodY = rev(method_order),
  stringsAsFactors = FALSE
) %>%
  mutate(
    ix = match(MethodX, method_order),
    iy = match(MethodY, method_order),
    value = case_when(
      ix < iy ~ spearman_mat[cbind(ix, iy)],
      ix > iy ~ jaccard_mat[cbind(ix, iy)],
      TRUE ~ 1
    ),
    label = sprintf("%.2f", value),
    MethodX = factor(MethodX, levels = method_order),
    MethodY = factor(MethodY, levels = rev(method_order))
  )

panel_c <- ggplot(matrix_long, aes(x = MethodX, y = MethodY, fill = value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = label), size = 3.1, color = "black", fontface = "bold") +
  scale_fill_gradient2(
    low = "#4A78A8",
    mid = "#F3F0F0",
    high = "#B30000",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Value"
  ) +
  labs(
    title = NULL,
    x = NULL,
    y = NULL,
    caption = "Upper triangle: Spearman correlation of -log10(p) | Lower triangle: Jaccard similarity of Bonferroni-significant sets"
  ) +
  panel_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.caption = element_text(size = 9, color = "#444444")
  )

# ------------------------------------------------------------------------------
# Panel D: Integrated candidate-gene prioritization scatter
# ------------------------------------------------------------------------------
gene_df <- gene_df %>%
  mutate(
    gwas_min_p = suppressWarnings(as.numeric(gwas_min_p)),
    magma_p = suppressWarnings(as.numeric(magma_p)),
    hit_gwas = as.logical(hit_gwas),
    hit_magma = as.logical(hit_magma),
    hit_pathway = as.logical(hit_pathway)
  ) %>%
  filter(is.finite(gwas_min_p), is.finite(magma_p), gwas_min_p > 0, magma_p > 0) %>%
  mutate(
    evidence = case_when(
      hit_gwas & hit_magma & hit_pathway ~ "All 3",
      hit_gwas & hit_pathway & !hit_magma ~ "GWAS + Pathway",
      hit_magma & hit_pathway & !hit_gwas ~ "MAGMA + Pathway",
      hit_gwas & hit_magma & !hit_pathway ~ "GWAS + MAGMA",
      hit_gwas & !hit_magma & !hit_pathway ~ "GWAS only",
      hit_magma & !hit_gwas & !hit_pathway ~ "MAGMA only",
      hit_pathway & !hit_gwas & !hit_magma ~ "Pathway only",
      TRUE ~ "None"
    ),
    evidence = factor(
      evidence,
      levels = c("All 3", "GWAS + Pathway", "MAGMA + Pathway", "GWAS + MAGMA",
                 "GWAS only", "MAGMA only", "Pathway only", "None")
    ),
    log_gwas = -log10(gwas_min_p),
    log_magma = -log10(magma_p),
    point_size = ifelse(evidence == "None", 0.85, 1.8)
  ) %>%
  arrange(evidence == "None", log_gwas, log_magma)

gwas_thresh <- max(gene_df$gwas_min_p[gene_df$hit_gwas], na.rm = TRUE)
magma_thresh <- max(gene_df$magma_p[gene_df$hit_magma], na.rm = TRUE)

panel_d <- ggplot(gene_df, aes(x = log_gwas, y = log_magma, color = evidence, size = point_size)) +
  geom_vline(xintercept = -log10(gwas_thresh), linetype = "dashed", color = "gray55", linewidth = 0.45) +
  geom_hline(yintercept = -log10(magma_thresh), linetype = "dashed", color = "gray55", linewidth = 0.45) +
  geom_point(alpha = 0.7) +
  scale_size_identity() +
  scale_color_manual(values = c(
    "All 3" = "#FF2B2B",
    "GWAS + Pathway" = "#FF66B3",
    "MAGMA + Pathway" = "#138D90",
    "GWAS + MAGMA" = "#A64DFF",
    "GWAS only" = "#F0A11E",
    "MAGMA only" = "#243B9B",
    "Pathway only" = "#3FA34D",
    "None" = "#D0D0D0"
  ), drop = FALSE) +
  labs(
    title = NULL,
    x = expression(-log[10](GWAS ~ italic(p))),
    y = expression(-log[10](MAGMA ~ italic(p)))
  ) +
  panel_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 9)
  ) +
  guides(
    color = guide_legend(
      nrow = 2,
      byrow = TRUE,
      override.aes = list(size = 3.5, alpha = 1)
    )
  )

# ------------------------------------------------------------------------------
# Assemble and save
# ------------------------------------------------------------------------------
fig7 <- (
  wrap_elements(full = panel_a_inner) |
    panel_b
) / (
  panel_c |
    panel_d
) +
  plot_layout(widths = c(1.05, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(face = "bold", size = 22))

ggsave(
  filename = file.path(final_fig_dir, "Fig7.png"),
  plot = fig7,
  width = 16,
  height = 14,
  dpi = 300,
  bg = "white"
)

message("Fig7 output written to Final/final_main_fig/Fig7.png")
