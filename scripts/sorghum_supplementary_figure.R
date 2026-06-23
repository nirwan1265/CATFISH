# ==============================================================================
# Supplementary Figure: CATFISH Component Test Comparison
# Sorghum Stem Volume
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(UpSetR)
library(cowplot)

OUT_DIR <- "sorghum_catfish_results"

# ==============================================================================
# CONFIGURATION - Change tau option here
# ==============================================================================
# Options: "default", "strict"
TAU_OPTION <- "strict"

# Set input file and output suffix based on tau option
if (TAU_OPTION == "strict") {
  PATHWAY_FILE <- file.path(OUT_DIR, "sorghum_stem_vol_CATFISH_B1000000_GPD_strict_tau.csv")
  OUT_SUFFIX <- "_strict_tau"
} else {
  PATHWAY_FILE <- file.path(OUT_DIR, "sorghum_stem_vol_CATFISH_B1000000_GPD.csv")
  OUT_SUFFIX <- ""
}
# ==============================================================================

# ==============================================================================
# Plot Theme (matching all_figs.R)
# ==============================================================================

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 20, color = "black"),
    axis.text.y = element_text(size = 20, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 20)
  )

# ==============================================================================
# Load CATFISH Results
# ==============================================================================

cat("=== Loading CATFISH Results ===\n")
cat("Using:", PATHWAY_FILE, "\n")

omni_results <- read.csv(PATHWAY_FILE, stringsAsFactors = FALSE)

cat("Total pathways:", nrow(omni_results), "\n")

# Define p-value columns
comp_cols <- c(
  "ACAT" = "acat_p",
  "Fisher" = "fisher_p",
  "TFisher" = "tfisher_p_analytic",
  "minP" = "minp_p_analytic",
  "Stouffer" = "stouffer_p_analytic"
)

omni_col <- "omni_p_final"

# ==============================================================================
# Panel A: UpSet Plot - Component Tests vs Omnibus
# ==============================================================================

cat("\n=== Creating Panel A: UpSet Plot ===\n")

# Bonferroni threshold
n_pathways <- nrow(omni_results)
bonf_threshold <- 0.05 / n_pathways

cat("Bonferroni threshold:", format(bonf_threshold, digits = 3), "\n")

# Get significant pathways for each method
get_sig <- function(col) {
  omni_results$pathway_id[omni_results[[col]] < bonf_threshold]
}

sig_sets <- list(
  "Omnibus" = get_sig(omni_col),
  "ACAT" = get_sig(comp_cols["ACAT"]),
  "Fisher" = get_sig(comp_cols["Fisher"]),
  "TFisher" = get_sig(comp_cols["TFisher"]),
  "minP" = get_sig(comp_cols["minP"]),
  "Stouffer" = get_sig(comp_cols["Stouffer"])
)

cat("Significant pathways per method:\n")
for (m in names(sig_sets)) {
  cat("  ", m, ":", length(sig_sets[[m]]), "\n")
}

# Create UpSet matrix
all_pathways <- unique(unlist(sig_sets))
upset_matrix <- sapply(names(sig_sets), function(s) {
  as.integer(all_pathways %in% sig_sets[[s]])
})
rownames(upset_matrix) <- all_pathways
upset_df <- as.data.frame(upset_matrix)

# Save UpSet plot
png(file.path(OUT_DIR, paste0("SuppFig_A_UpSet", OUT_SUFFIX, ".png")),
    width = 3600, height = 2800, res = 300)
upset(
  upset_df,
  sets = c("Omnibus", "ACAT", "Fisher", "TFisher", "minP", "Stouffer"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(2.5, 2.5, 2.5, 2.5, 2.5, 2.5),
  point.size = 6,
  line.size = 2,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = c("darkred", rep("steelblue", 5)),
  main.bar.color = "gray30",
  set_size.angles = 45
)
dev.off()

cat("Panel A saved\n")

# ==============================================================================
# Panel B: QQ Plot - Omnibus P-values with BH Significance
# ==============================================================================

cat("\n=== Creating Panel B: Omnibus QQ Plot ===\n")

# Get p-values for all methods (for later use)
all_methods <- c(comp_cols, "Omnibus" = omni_col)

# Use omnibus p-values for QQ plot
pvals_omni <- omni_results$omni_p_final
pvals_omni <- pvals_omni[!is.na(pvals_omni) & pvals_omni > 0 & pvals_omni < 1]
pvals_sorted <- sort(pvals_omni)
n_pw <- length(pvals_sorted)

# Expected under null
expected_p <- (1:n_pw) / (n_pw + 1)

qq_df <- data.frame(
  expected = -log10(expected_p),
  observed = -log10(pvals_sorted)
)

# Calculate thresholds
bh_pvals <- p.adjust(sort(omni_results$omni_p_final), method = "BH")
bh_05_idx <- max(which(bh_pvals < 0.05))
bh_05_pval <- sort(omni_results$omni_p_final)[bh_05_idx]

# Color points by significance
qq_df$significance <- case_when(
  pvals_sorted < bonf_threshold ~ "Bonferroni < 0.05",
  pvals_sorted < bh_05_pval ~ "BH < 0.05",
  TRUE ~ "Not significant"
)
qq_df$significance <- factor(qq_df$significance,
                             levels = c("Bonferroni < 0.05", "BH < 0.05", "Not significant"))

n_bonf <- sum(pvals_sorted < bonf_threshold)
n_bh <- sum(bh_pvals < 0.05)

cat("Pathways Bonferroni < 0.05:", n_bonf, "\n")
cat("Pathways BH < 0.05:", n_bh, "\n")

panel_b <- ggplot(qq_df, aes(x = expected, y = observed, color = significance)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
  geom_point(alpha = 0.7, size = 2.5) +
  scale_color_manual(values = c(
    "Bonferroni < 0.05" = "darkred",
    "BH < 0.05" = "steelblue",
    "Not significant" = "gray60"
  )) +
  labs(
    x = expression(Expected ~ -log[10](p)),
    y = expression(Observed ~ -log[10](p)),
    caption = paste0("n = ", n_pw, " pathways | Bonferroni: ", n_bonf, " | BH < 0.05: ", n_bh)
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 16),
    legend.title = element_blank(),
    plot.caption = element_text(size = 14, hjust = 0.5)
  ) +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(file.path(OUT_DIR, paste0("SuppFig_B_QQplot", OUT_SUFFIX, ".png")), panel_b,
       width = 10, height = 9, dpi = 300, bg = "white")

cat("Panel B saved\n")

# ==============================================================================
# Panel C: Pairwise Jaccard Similarity and Correlation
# ==============================================================================

cat("\n=== Creating Panel C: Jaccard and Correlation Matrix ===\n")

method_names <- c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")
n_methods <- length(method_names)

# Calculate Jaccard similarity between significant pathway sets
jaccard_mat <- matrix(NA, n_methods, n_methods,
                      dimnames = list(method_names, method_names))

for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    set_i <- sig_sets[[method_names[i]]]
    set_j <- sig_sets[[method_names[j]]]
    inter <- length(intersect(set_i, set_j))
    uni <- length(union(set_i, set_j))
    jaccard_mat[i, j] <- if (uni > 0) inter / uni else 0
  }
}

# Calculate correlation between p-values
cor_mat <- matrix(NA, n_methods, n_methods,
                  dimnames = list(method_names, method_names))

for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    col_i <- all_methods[method_names[i]]
    col_j <- all_methods[method_names[j]]
    pvals_i <- -log10(omni_results[[col_i]])
    pvals_j <- -log10(omni_results[[col_j]])
    valid <- !is.na(pvals_i) & !is.na(pvals_j) & is.finite(pvals_i) & is.finite(pvals_j)
    cor_mat[i, j] <- cor(pvals_i[valid], pvals_j[valid], method = "spearman")
  }
}

# Create combined matrix: upper triangle = correlation, lower triangle = Jaccard
combined_mat <- matrix(NA, n_methods, n_methods,
                       dimnames = list(method_names, method_names))
for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    if (i < j) {
      combined_mat[i, j] <- cor_mat[i, j]  # Upper: correlation
    } else if (i > j) {
      combined_mat[i, j] <- jaccard_mat[i, j]  # Lower: Jaccard
    } else {
      combined_mat[i, j] <- 1  # Diagonal
    }
  }
}

# Convert to long format for plotting
combined_long <- as.data.frame(as.table(combined_mat))
names(combined_long) <- c("Method1", "Method2", "Value")
combined_long$Method1 <- factor(combined_long$Method1, levels = method_names)
combined_long$Method2 <- factor(combined_long$Method2, levels = method_names)

# Add type indicator
combined_long$Type <- ifelse(
  as.numeric(combined_long$Method1) < as.numeric(combined_long$Method2),
  "Correlation",
  ifelse(as.numeric(combined_long$Method1) > as.numeric(combined_long$Method2),
         "Jaccard", "Diagonal")
)

# Create labels
combined_long$Label <- sprintf("%.2f", combined_long$Value)
combined_long$Label[combined_long$Type == "Diagonal"] <- ""

panel_c <- ggplot(combined_long, aes(x = Method2, y = Method1, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 5, fontface = "bold") +
  scale_fill_gradient2(
    low = "white", mid = "steelblue", high = "darkred",
    midpoint = 0.5, limits = c(0, 1),
    name = "Value"
  ) +
  labs(
    x = "", y = "",
    caption = "Upper triangle: Spearman correlation of -log10(p)\nLower triangle: Jaccard similarity of significant sets"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    legend.position = "right",
    legend.title = element_text(size = 16, face = "bold"),
    plot.caption = element_text(size = 14, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed()

ggsave(file.path(OUT_DIR, paste0("SuppFig_C_Jaccard_Corr", OUT_SUFFIX, ".png")), panel_c,
       width = 10, height = 9, dpi = 300, bg = "white")

cat("Panel C saved\n")

# ==============================================================================
# Panel D: P-value Distribution (Box Plot)
# ==============================================================================

cat("\n=== Creating Panel D: P-value Distribution ===\n")

# Prepare data for box plot
pval_long <- omni_results %>%
  select(pathway_id, all_of(unname(all_methods))) %>%
  pivot_longer(
    cols = -pathway_id,
    names_to = "method_col",
    values_to = "p_value"
  ) %>%
  mutate(
    method = case_when(
      method_col == "acat_p" ~ "ACAT",
      method_col == "fisher_p" ~ "Fisher",
      method_col == "tfisher_p_analytic" ~ "TFisher",
      method_col == "minp_p_analytic" ~ "minP",
      method_col == "stouffer_p_analytic" ~ "Stouffer",
      method_col == "omni_p_final" ~ "Omnibus",
      TRUE ~ method_col
    ),
    method = factor(method, levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus")),
    neg_log10_p = -log10(pmax(p_value, 1e-16))
  )

panel_d <- ggplot(pval_long, aes(x = method, y = neg_log10_p, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 1) +
  geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed",
             color = "red", linewidth = 0.8) +
  scale_fill_manual(values = c(
    "ACAT" = "#E41A1C", "Fisher" = "#377EB8", "TFisher" = "#4DAF4A",
    "minP" = "#984EA3", "Stouffer" = "#FF7F00", "Omnibus" = "#000000"
  )) +
  labs(
    x = "Method",
    y = expression(-log[10](p)),
    caption = paste0("Red dashed line: Bonferroni threshold (-log10(",
                     format(bonf_threshold, digits = 2), ") = ",
                     round(-log10(bonf_threshold), 2), ")")
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 20),
    legend.position = "none",
    plot.caption = element_text(size = 14, hjust = 0.5)
  )

ggsave(file.path(OUT_DIR, paste0("SuppFig_D_Pval_Dist", OUT_SUFFIX, ".png")), panel_d,
       width = 10, height = 8, dpi = 300, bg = "white")

cat("Panel D saved\n")

# ==============================================================================
# Combine All Panels
# ==============================================================================

cat("\n=== Combining All Panels ===\n")

# Read Panel A as image
panel_a_plot <- ggdraw() +
  draw_image(file.path(OUT_DIR, paste0("SuppFig_A_UpSet", OUT_SUFFIX, ".png"))) +
  draw_label("A", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

panel_b_plot <- ggdraw(panel_b) +
  draw_label("B", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

panel_c_plot <- ggdraw(panel_c) +
  draw_label("C", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

panel_d_plot <- ggdraw(panel_d) +
  draw_label("D", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

# Combine using cowplot
top_row <- plot_grid(panel_a_plot, panel_b_plot, ncol = 2,
                     rel_widths = c(1.2, 1), align = "h")
bottom_row <- plot_grid(panel_c_plot, panel_d_plot, ncol = 2,
                        rel_widths = c(1, 1), align = "h")

supp_fig <- plot_grid(top_row, bottom_row, nrow = 2,
                      rel_heights = c(1.1, 1))

ggsave(file.path(OUT_DIR, paste0("Supplementary_Figure_Component_Tests", OUT_SUFFIX, ".png")), supp_fig,
       width = 22, height = 20, dpi = 300, bg = "white")

ggsave(file.path(OUT_DIR, paste0("Supplementary_Figure_Component_Tests", OUT_SUFFIX, ".pdf")), supp_fig,
       width = 22, height = 20, bg = "white")

cat("\n=== DONE ===\n")
cat("Supplementary Figure saved:\n")
cat("  -", file.path(OUT_DIR, "Supplementary_Figure_Component_Tests.png"), "\n")
cat("  -", file.path(OUT_DIR, "Supplementary_Figure_Component_Tests.pdf"), "\n")
cat("\nIndividual panels:\n")
cat("  -", file.path(OUT_DIR, "SuppFig_A_UpSet.png"), "\n")
cat("  -", file.path(OUT_DIR, "SuppFig_B_QQplot.png"), "\n")
cat("  -", file.path(OUT_DIR, "SuppFig_C_Jaccard_Corr.png"), "\n")
cat("  -", file.path(OUT_DIR, "SuppFig_D_Pval_Dist.png"), "\n")

# Print summary stats
cat("\n=== Summary Statistics ===\n")
cat("\nJaccard similarity matrix:\n")
print(round(jaccard_mat, 3))
cat("\nCorrelation matrix:\n")
print(round(cor_mat, 3))
