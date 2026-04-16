# ==============================================================================
# Compare Component Tests from CATFISH Omnibus Results - Sorghum Stem Volume
# Includes QQ plot with broken axis for B=1M GPD results
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# For Venn diagrams
if (!requireNamespace("ggVennDiagram", quietly = TRUE)) {
  install.packages("ggVennDiagram")
}
library(ggVennDiagram)

# ==============================================================================
# Output directory
# ==============================================================================
fig_dir <- "Figures/Fig_Sorghum_Component_Tests"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

ind_fig_dir <- file.path(fig_dir, "Ind_figs")
if (!dir.exists(ind_fig_dir)) dir.create(ind_fig_dir, recursive = TRUE)

# ==============================================================================
# Load Results
# ==============================================================================

omni_results <- read.csv("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv",
                         stringsAsFactors = FALSE)

n_pathways <- nrow(omni_results)
bonf_thresh <- 0.05 / n_pathways

cat("Loaded", n_pathways, "pathways\n")
cat("Bonferroni threshold:", format(bonf_thresh, digits = 3), "\n")

# ==============================================================================
# Define Component P-value Columns
# ==============================================================================

# Use analytic component p-values
use_cols <- c(
  "ACAT"     = "acat_p",
  "Fisher"   = "fisher_p",
  "TFisher"  = "tfisher_p_analytic",
  "minP"     = "minp_p_analytic",
  "Stouffer" = "stouffer_p_analytic"
)

# Add omnibus
use_cols <- c(use_cols, "Omnibus" = "omni_p_final")

# Check counts at Bonferroni
cat("\nSignificant pathways at Bonferroni level:\n")
for (method in names(use_cols)) {
  col <- use_cols[method]
  if (col %in% names(omni_results)) {
    n_sig <- sum(omni_results[[col]] < bonf_thresh, na.rm = TRUE)
    cat(sprintf("  %s: %d\n", method, n_sig))
  }
}

# ==============================================================================
# Create Significant Pathway Sets
# ==============================================================================

get_significant_pathways <- function(df, col, threshold) {
  if (!col %in% names(df)) return(character(0))
  pvals <- df[[col]]
  pathway_ids <- df$pathway_id
  valid_idx <- !is.na(pvals)
  pathway_ids[valid_idx][pvals[valid_idx] < threshold]
}

sig_pathways <- lapply(names(use_cols), function(method) {
  get_significant_pathways(omni_results, use_cols[method], bonf_thresh)
})
names(sig_pathways) <- names(use_cols)

# ==============================================================================
# PANEL A: QQ Plot with Broken Axis
# ==============================================================================

qq_data <- omni_results %>%
  arrange(omni_p_final) %>%
  mutate(
    rank = row_number(),
    expected_p = rank / (n_pathways + 1),
    expected_neglog10 = -log10(expected_p),
    observed_neglog10 = -log10(omni_p_final),
    significance = ifelse(omni_p_final < bonf_thresh, "Bonferroni significant", "Non-significant")
  )

qq_data$significance <- factor(qq_data$significance, levels = c("Bonferroni significant", "Non-significant"))

n_sig <- sum(omni_results$omni_p_final < bonf_thresh)
max_obs <- max(qq_data$observed_neglog10)

# Transform y-axis: compress values > 5
qq_broken <- qq_data %>%
  mutate(
    y_display = case_when(
      observed_neglog10 <= 5 ~ observed_neglog10,
      observed_neglog10 > 5 ~ 5 + (observed_neglog10 - 5) * 0.12
    )
  )

max_y_display <- max(qq_broken$y_display)
y_breaks <- c(0, 1, 2, 3, 4, 5, 5 + (10-5)*0.12, 5 + (20-5)*0.12, 5 + (30-5)*0.12)
y_labels <- c("0", "1", "2", "3", "4", "5", "10", "20", "30")

panelA <- ggplot(qq_broken, aes(x = expected_neglog10, y = y_display)) +
  annotate("segment", x = 0, xend = 3, y = 0, yend = 3, linetype = "dashed", color = "gray50") +
  annotate("rect", xmin = -0.15, xmax = 3.05, ymin = 4.95, ymax = 5.15, fill = "white", color = NA) +
  annotate("segment", x = -0.05, xend = 0.1, y = 4.9, yend = 5.2, color = "gray30", linewidth = 0.5) +
  annotate("segment", x = 0.1, xend = 0.25, y = 4.9, yend = 5.2, color = "gray30", linewidth = 0.5) +
  geom_point(aes(color = significance), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Bonferroni significant" = "#D7191C", "Non-significant" = "gray60"), name = NULL) +
  geom_hline(yintercept = -log10(bonf_thresh), linetype = "dotted", color = "#D7191C", linewidth = 0.6) +
  scale_y_continuous(breaks = y_breaks, labels = y_labels) +
  labs(x = expression(Expected ~ -log[10](P)), y = expression(Observed ~ -log[10](P)), tag = "A") +
  coord_cartesian(xlim = c(0, 3), ylim = c(0, max_y_display + 0.2), clip = "off") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = c(0.7, 0.2),
    legend.background = element_rect(fill = "white", color = "gray80"),
    legend.text = element_text(size = 8),
    panel.grid.minor = element_blank(),
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "qq_omnibus_broken.png"), panelA, width = 5, height = 5, dpi = 300)

# ==============================================================================
# PANEL B: Venn Diagram (5 component tests)
# ==============================================================================

# Use only component tests (not omnibus) for Venn
sig_components <- sig_pathways[names(sig_pathways) != "Omnibus"]

venn_plot <- ggVennDiagram(sig_components, label = "count", label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "steelblue", name = "Count") +
  labs(tag = "B") +
  theme_void(base_size = 11) +
  theme(
    legend.position = "none",
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "venn_component_tests.png"), venn_plot, width = 6, height = 5, dpi = 300)

# ==============================================================================
# PANEL C: Jaccard Similarity Heatmap
# ==============================================================================

jaccard <- function(a, b) {
  if (length(a) == 0 && length(b) == 0) return(1)
  inter <- length(intersect(a, b))
  union_size <- length(union(a, b))
  if (union_size == 0) return(0)
  inter / union_size
}

methods <- names(sig_components)
n_methods <- length(methods)

jaccard_matrix <- matrix(NA, nrow = n_methods, ncol = n_methods, dimnames = list(methods, methods))
for (i in 1:n_methods) {
  for (j in 1:n_methods) {
    jaccard_matrix[i, j] <- jaccard(sig_components[[i]], sig_components[[j]])
  }
}

jaccard_df <- as.data.frame(as.table(jaccard_matrix))
names(jaccard_df) <- c("Method1", "Method2", "Jaccard")

panelC <- ggplot(jaccard_df, aes(x = Method1, y = Method2, fill = Jaccard)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Jaccard)), color = "black", size = 3) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue", midpoint = 0.5, limits = c(0, 1)) +
  labs(x = NULL, y = NULL, tag = "C") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "jaccard_heatmap.png"), panelC, width = 5, height = 4, dpi = 300)

# ==============================================================================
# PANEL D: Spearman Rank Correlation Heatmap
# ==============================================================================

rank_df <- data.frame(pathway_id = omni_results$pathway_id)
for (method in names(use_cols)) {
  col <- use_cols[method]
  if (col %in% names(omni_results)) {
    rank_df[[method]] <- rank(omni_results[[col]], ties.method = "average", na.last = "keep")
  }
}

rank_matrix <- as.matrix(rank_df[, -1])
spearman_cor <- cor(rank_matrix, use = "pairwise.complete.obs", method = "spearman")

spearman_df <- as.data.frame(as.table(spearman_cor))
names(spearman_df) <- c("Method1", "Method2", "Correlation")

panelD <- ggplot(spearman_df, aes(x = Method1, y = Method2, fill = Correlation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.2f", Correlation)), color = "black", size = 3) +
  scale_fill_gradient2(low = "white", mid = "steelblue", high = "darkblue", midpoint = 0.5, limits = c(0, 1)) +
  labs(x = NULL, y = NULL, tag = "D") +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "spearman_heatmap.png"), panelD, width = 6, height = 5, dpi = 300)

# ==============================================================================
# PANEL E: P-value Boxplot
# ==============================================================================

pval_long <- omni_results %>%
  select(pathway_id, all_of(unname(use_cols[use_cols %in% names(omni_results)]))) %>%
  pivot_longer(cols = -pathway_id, names_to = "col", values_to = "pvalue") %>%
  mutate(
    method = case_when(
      col == "acat_p" ~ "ACAT",
      col == "fisher_p" ~ "Fisher",
      col == "tfisher_p_analytic" ~ "TFisher",
      col == "minp_p_analytic" ~ "minP",
      col == "stouffer_p_analytic" ~ "Stouffer",
      col == "omni_p_final" ~ "Omnibus",
      TRUE ~ col
    ),
    log10_pval = -log10(pvalue)
  )

pval_long$method <- factor(pval_long$method, levels = c("ACAT", "Fisher", "TFisher", "minP", "Stouffer", "Omnibus"))

panelE <- ggplot(pval_long, aes(x = method, y = log10_pval, fill = method)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  geom_hline(yintercept = -log10(bonf_thresh), linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  labs(x = NULL, y = expression(-log[10](P)), tag = "E") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "pvalue_boxplot.png"), panelE, width = 5, height = 4, dpi = 300)

# ==============================================================================
# PANEL F: P-value Density
# ==============================================================================

panelF <- ggplot(pval_long, aes(x = log10_pval, fill = method, color = method)) +
  geom_density(alpha = 0.3) +
  geom_vline(xintercept = -log10(bonf_thresh), linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(x = expression(-log[10](P)), y = "Density", tag = "F", fill = "Method", color = "Method") +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),
    plot.tag = element_text(size = 14, face = "bold")
  )

ggsave(file.path(ind_fig_dir, "pvalue_density.png"), panelF, width = 6, height = 4, dpi = 300)

# ==============================================================================
# Combine All Panels
# ==============================================================================

# Layout:
# Row 1: A (QQ) | B (Venn)
# Row 2: C (Jaccard) | D (Spearman)
# Row 3: E (Boxplot) | F (Density)

combined_fig <- (panelA | venn_plot) / (panelC | panelD) / (panelE | panelF) +
  plot_layout(heights = c(1, 1, 0.8))

ggsave(file.path(fig_dir, "Fig_Component_Tests_Sorghum.png"), combined_fig, width = 12, height = 14, dpi = 300)
ggsave(file.path(fig_dir, "Fig_Component_Tests_Sorghum.pdf"), combined_fig, width = 12, height = 14)

cat("\n========================================\n")
cat("FIGURES SAVED\n")
cat("========================================\n")
cat("Combined figure:", file.path(fig_dir, "Fig_Component_Tests_Sorghum.png/pdf"), "\n")
cat("Individual panels:", ind_fig_dir, "\n")

# ==============================================================================
# Summary Statistics
# ==============================================================================

cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
cat("Total pathways:", n_pathways, "\n")
cat("Bonferroni threshold:", format(bonf_thresh, digits = 3), "\n")
cat("Omnibus significant:", n_sig, "\n")
cat("Max -log10(P):", round(max_obs, 1), "\n")
cat("\nAverage Jaccard similarity:", round(mean(jaccard_matrix[lower.tri(jaccard_matrix)]), 3), "\n")
cat("Average Spearman correlation:", round(mean(spearman_cor[lower.tri(spearman_cor)]), 3), "\n")
