# Pathway and Gene Comparison: CATFISH Omnibus vs MAGMA
# Compare top pathways and gene scores between methods

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# =============================================================================
# Load Results
# =============================================================================

# CATFISH omnibus results (B=1M GPD)
catfish <- read.csv("sorghum_catfish_results/sorghum_stem_vol_CATFISH_B1000000_GPD.csv",
                    stringsAsFactors = FALSE)

# MAGMA pathway results
magma_raw <- read.table("sorghum_catfish_results/sorghum_stem_vol_magma_pathway.gsa.out",
                        header = TRUE, comment.char = "#", stringsAsFactors = FALSE)

# Scored genes (CATFISH scoring)
scored_genes <- read.csv("sorghum_catfish_results/candidate_genes_GPD_B1M_scored.csv",
                         stringsAsFactors = FALSE)

# Gene-level MAGMA results
gene_results <- read.table("sorghum_catfish_results/sorghum_stem_vol_genes_combined.txt",
                           header = TRUE, stringsAsFactors = FALSE)

# =============================================================================
# Prepare data
# =============================================================================

# Clean up MAGMA pathway results
magma <- magma_raw %>%
  select(pathway_id = VARIABLE, magma_p = P, magma_ngenes = NGENES, magma_beta = BETA) %>%
  mutate(magma_neglog10p = -log10(magma_p))

# Clean up CATFISH results
catfish_clean <- catfish %>%
  select(pathway_id, pathway_name, n_genes, omni_p_final) %>%
  mutate(catfish_neglog10p = -log10(omni_p_final))

# Merge the two
comparison <- merge(catfish_clean, magma, by = "pathway_id", all.x = TRUE)

# Bonferroni threshold
bonf_thresh <- 0.05 / nrow(comparison)

# Add significance flags
comparison <- comparison %>%
  mutate(
    sig_catfish = omni_p_final < bonf_thresh,
    sig_magma = magma_p < bonf_thresh
  )

cat("=============================================================================\n")
cat("PATHWAY COMPARISON: CATFISH vs MAGMA\n")
cat("=============================================================================\n\n")

# =============================================================================
# PANEL A: Top Pathways with CATFISH-only, Common, MAGMA-only blocks
# =============================================================================

# Get top 10 from each method
top10_catfish_ids <- comparison %>% arrange(omni_p_final) %>% head(10) %>% pull(pathway_id)
top10_magma_ids <- comparison %>% arrange(magma_p) %>% head(10) %>% pull(pathway_id)

# Categorize
common_ids <- intersect(top10_catfish_ids, top10_magma_ids)
catfish_only_ids <- setdiff(top10_catfish_ids, top10_magma_ids)
magma_only_ids <- setdiff(top10_magma_ids, top10_catfish_ids)

cat(sprintf("Top 10 overlap: %d pathways\n", length(common_ids)))
cat(sprintf("CATFISH-only in top 10: %d pathways\n", length(catfish_only_ids)))
cat(sprintf("MAGMA-only in top 10: %d pathways\n", length(magma_only_ids)))

# Build plotting data
plot_pathways <- comparison %>%
  filter(pathway_id %in% c(catfish_only_ids, common_ids, magma_only_ids)) %>%
  mutate(
    category = case_when(
      pathway_id %in% common_ids ~ "Common",
      pathway_id %in% catfish_only_ids ~ "CATFISH Omnibus",
      pathway_id %in% magma_only_ids ~ "MAGMA"
    ),
    category = factor(category, levels = c("CATFISH Omnibus", "Common", "MAGMA"))
  )

# Create short names
plot_pathways$short_name <- sapply(plot_pathways$pathway_name, function(x) {
  x <- gsub("<[^>]+>", "", x)  # Remove HTML tags
  if (nchar(x) > 35) paste0(substr(x, 1, 32), "...") else x
})

# Order: CATFISH-only by catfish p, then common, then MAGMA-only by magma p
catfish_order <- plot_pathways %>%
  filter(category == "CATFISH Omnibus") %>%
  arrange(omni_p_final) %>%
  pull(pathway_id)
common_order <- plot_pathways %>%
  filter(category == "Common") %>%
  arrange(omni_p_final) %>%
  pull(pathway_id)
magma_order <- plot_pathways %>%
  filter(category == "MAGMA") %>%
  arrange(magma_p) %>%
  pull(pathway_id)

pathway_order <- c(catfish_order, common_order, magma_order)
plot_pathways$pathway_id <- factor(plot_pathways$pathway_id, levels = rev(pathway_order))

# Prepare long format for side-by-side bars
plot_long <- plot_pathways %>%
  select(pathway_id, short_name, n_genes, category, catfish_neglog10p, magma_neglog10p) %>%
  pivot_longer(cols = c(catfish_neglog10p, magma_neglog10p),
               names_to = "method", values_to = "neglog10p") %>%
  mutate(method = ifelse(method == "catfish_neglog10p", "CATFISH", "MAGMA"))

# Panel A - using facets for vertical block labels
plot_long$category <- factor(plot_long$category, levels = c("CATFISH Omnibus", "Common", "MAGMA"))

# Reorder within each category
plot_long <- plot_long %>%
  group_by(category) %>%
  arrange(category, desc(neglog10p)) %>%
  ungroup()

panelA <- ggplot(plot_long, aes(x = neglog10p, y = pathway_id, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_vline(xintercept = -log10(bonf_thresh), linetype = "dashed", color = "red", alpha = 0.7) +
  scale_fill_manual(values = c("CATFISH" = "#2166AC", "MAGMA" = "#B2182B")) +
  scale_y_discrete(labels = function(x) {
    sapply(x, function(id) {
      row <- plot_pathways[plot_pathways$pathway_id == id, ]
      if (nrow(row) > 0) {
        sprintf("%s (n=%d)", row$short_name[1], row$n_genes[1])
      } else id
    })
  }) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = expression(-log[10](p-value)),
    y = NULL,
    fill = "Method",
    tag = "A"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 9),
    strip.background = element_rect(fill = "gray90", color = "gray50"),
    panel.spacing = unit(0.3, "lines"),
    plot.tag = element_text(size = 14, face = "bold")
  )

# =============================================================================
# PANEL B: Top 10 Genes - CATFISH Score vs MAGMA
# =============================================================================

# Top 10 by CATFISH score
top10_catfish_genes <- scored_genes %>%
  arrange(desc(score)) %>%
  head(10) %>%
  pull(GENE)

# Top 10 by MAGMA gene p-value
top10_magma_genes <- gene_results %>%
  arrange(P) %>%
  head(10) %>%
  pull(GENE)

# Categorize genes
common_genes <- intersect(top10_catfish_genes, top10_magma_genes)
catfish_only_genes <- setdiff(top10_catfish_genes, top10_magma_genes)
magma_only_genes <- setdiff(top10_magma_genes, top10_catfish_genes)

cat(sprintf("\nTop 10 genes overlap: %d genes\n", length(common_genes)))

# Build gene plotting data
gene_data <- scored_genes %>%
  filter(GENE %in% c(catfish_only_genes, common_genes, magma_only_genes)) %>%
  mutate(
    category = case_when(
      GENE %in% common_genes ~ "Common",
      GENE %in% catfish_only_genes ~ "CATFISH Score",
      GENE %in% magma_only_genes ~ "MAGMA"
    ),
    category = factor(category, levels = c("CATFISH Score", "Common", "MAGMA")),
    # Create score for both methods (normalized)
    catfish_score = score,
    magma_score = -log10(magma_p)  # Use -log10(p) as MAGMA "score"
  )

# Order genes
catfish_gene_order <- gene_data %>%
  filter(category == "CATFISH Score") %>%
  arrange(desc(score)) %>%
  pull(GENE)
common_gene_order <- gene_data %>%
  filter(category == "Common") %>%
  arrange(desc(score)) %>%
  pull(GENE)
magma_gene_order <- gene_data %>%
  filter(category == "MAGMA") %>%
  arrange(magma_p) %>%
  pull(GENE)

gene_order <- c(catfish_gene_order, common_gene_order, magma_gene_order)
gene_data$GENE <- factor(gene_data$GENE, levels = rev(gene_order))

# Normalize scores to same scale for visualization
max_catfish <- max(gene_data$catfish_score, na.rm = TRUE)
max_magma <- max(gene_data$magma_score, na.rm = TRUE)

gene_long <- gene_data %>%
  mutate(
    catfish_norm = catfish_score / max_catfish * 10,
    magma_norm = magma_score / max_magma * 10
  ) %>%
  select(GENE, gene_category = category, catfish_norm, magma_norm) %>%
  pivot_longer(cols = c(catfish_norm, magma_norm),
               names_to = "method", values_to = "norm_score") %>%
  mutate(method = ifelse(method == "catfish_norm", "CATFISH Score", "MAGMA -log10(p)"))

# Panel B - using facets for vertical block labels
gene_long$gene_category <- factor(gene_long$gene_category, levels = c("CATFISH Score", "Common", "MAGMA"))

panelB <- ggplot(gene_long, aes(x = norm_score, y = GENE, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = c("CATFISH Score" = "#2166AC", "MAGMA -log10(p)" = "#B2182B")) +
  facet_grid(gene_category ~ ., scales = "free_y", space = "free_y") +
  labs(
    x = "Normalized Score (0-10 scale)",
    y = NULL,
    fill = "Method",
    tag = "B"
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "top",
    axis.text.y = element_text(size = 8),
    strip.text.y.right = element_text(angle = 0, face = "bold", size = 9),
    strip.background = element_rect(fill = "gray90", color = "gray50"),
    panel.spacing = unit(0.3, "lines"),
    plot.tag = element_text(size = 14, face = "bold")
  )

# =============================================================================
# PANEL C: Scatter plot of pathway rankings
# =============================================================================

comparison <- comparison %>%
  arrange(omni_p_final) %>%
  mutate(catfish_rank = row_number()) %>%
  arrange(magma_p) %>%
  mutate(magma_rank = row_number())

panelC <- ggplot(comparison, aes(x = catfish_rank, y = magma_rank)) +
  geom_point(aes(color = pmin(catfish_neglog10p, magma_neglog10p)), alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_viridis_c(name = expression(min(-log[10](p))), option = "C") +
  labs(
    x = "CATFISH Omnibus Rank",
    y = "MAGMA Rank",
    tag = "C"
  ) +
  annotate("text", x = 350, y = 50,
           label = sprintf("Spearman rho = %.3f",
                          cor(comparison$catfish_rank, comparison$magma_rank,
                              method = "spearman", use = "complete.obs")),
           size = 4, hjust = 1) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "right",
    plot.tag = element_text(size = 14, face = "bold")
  )

# =============================================================================
# Combine and Save Figures
# =============================================================================

# Create output directory
fig_dir <- "sorghum_catfish_results/figures"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# Combined figure
combined_fig <- (panelA / panelB) | panelC
combined_fig <- combined_fig + plot_layout(widths = c(2, 1))

ggsave(file.path(fig_dir, "pathway_omnibus_magma_comparison.png"),
       combined_fig, width = 14, height = 12, dpi = 300)
ggsave(file.path(fig_dir, "pathway_omnibus_magma_comparison.pdf"),
       combined_fig, width = 14, height = 12)

cat("\n\nFigures saved to:\n")
cat(sprintf("  %s/pathway_omnibus_magma_comparison.png\n", fig_dir))
cat(sprintf("  %s/pathway_omnibus_magma_comparison.pdf\n", fig_dir))

# Individual panels
ggsave(file.path(fig_dir, "panelA_pathway_comparison.png"), panelA, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "panelB_gene_comparison.png"), panelB, width = 10, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "panelC_rank_scatter.png"), panelC, width = 6, height = 6, dpi = 300)

# =============================================================================
# Summary Statistics
# =============================================================================

cat("\n\n=============================================================================\n")
cat("SUMMARY\n")
cat("=============================================================================\n")
cat(sprintf("Total pathways analyzed: %d\n", nrow(comparison)))
cat(sprintf("Bonferroni threshold: %.2e (0.05/%d)\n", bonf_thresh, nrow(comparison)))
cat(sprintf("\nPathway Results:\n"))
cat(sprintf("  Top 10 overlap: %d pathways\n", length(common_ids)))
cat(sprintf("  CATFISH-only in top 10: %d\n", length(catfish_only_ids)))
cat(sprintf("  MAGMA-only in top 10: %d\n", length(magma_only_ids)))
cat(sprintf("  Significant (CATFISH): %d\n", sum(comparison$sig_catfish)))
cat(sprintf("  Significant (MAGMA): %d\n", sum(comparison$sig_magma, na.rm = TRUE)))
cat(sprintf("\nGene Results:\n"))
cat(sprintf("  Top 10 overlap: %d genes\n", length(common_genes)))
cat(sprintf("  CATFISH-only in top 10: %d\n", length(catfish_only_genes)))
cat(sprintf("  MAGMA-only in top 10: %d\n", length(magma_only_genes)))
cat(sprintf("\nRank correlation (Spearman): %.4f\n",
            cor(comparison$catfish_rank, comparison$magma_rank, method = "spearman", use = "complete.obs")))

# Save comparison tables
write.csv(comparison, file.path(fig_dir, "pathway_comparison_table.csv"), row.names = FALSE)
write.csv(gene_data, file.path(fig_dir, "gene_comparison_table.csv"), row.names = FALSE)
cat(sprintf("\nTables saved to: %s/\n", fig_dir))
