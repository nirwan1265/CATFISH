# ==============================================================================
# Figure 5: Candidate Gene Analysis
# Sorghum Stem Volume
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)
library(patchwork)

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

if (!requireNamespace("cowplot", quietly = TRUE)) {
  install.packages("cowplot")
}
library(cowplot)

OUT_DIR <- "sorghum_catfish_results"

# ==============================================================================
# Plot Theme
# ==============================================================================

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x = element_text(size = 24, face = "bold"),
    axis.title.y = element_text(size = 24, face = "bold"),
    axis.text.x = element_text(size = 24, color = "black"),
    axis.text.y = element_text(size = 24, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 24)
  )

# ==============================================================================
# Load Data
# ==============================================================================

cat("=== Loading data ===\n")

# Load gene location
gene_loc <- read.delim("inst/extdata/sorghum.genes.loc", stringsAsFactors = FALSE)
gene_loc <- gene_loc %>%
  mutate(
    START_EXT = pmax(0, START - 25000),
    STOP_EXT = STOP + 25000
  )

# Load GWAS
gwas_raw <- fread(
  "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt",
  header = TRUE
)
names(gwas_raw)[names(gwas_raw) == "p_wald"] <- "P"
names(gwas_raw)[names(gwas_raw) == "ps"] <- "POS"
names(gwas_raw)[names(gwas_raw) == "chr"] <- "CHR"

# Map GWAS to genes
gwas_dt <- as.data.table(gwas_raw)
gene_dt <- as.data.table(gene_loc)
setkey(gene_dt, CHR, START_EXT, STOP_EXT)
gwas_dt[, c("START_EXT", "STOP_EXT") := .(POS, POS)]
setkey(gwas_dt, CHR, START_EXT, STOP_EXT)

snp_gene_map <- foverlaps(
  gwas_dt, gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within", nomatch = NULL
)

gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_logP = -log10(min(P, na.rm = TRUE)),
    gwas_n_snps = n(),
    .groups = "drop"
  ) %>%
  arrange(gwas_min_p)

# GWAS rank - same p-value gets same rank
gwas_gene <- gwas_gene %>%
  mutate(gwas_rank = dense_rank(gwas_min_p))

cat("GWAS genes:", nrow(gwas_gene), "\n")
cat("GWAS genes with -log10(p) > 7:", sum(gwas_gene$gwas_logP > 7), "\n")

# Load MAGMA
magma_gene <- read.table(
  file.path(OUT_DIR, "sorghum_stem_vol_genes_combined.txt"),
  header = TRUE, stringsAsFactors = FALSE
)
magma_gene <- magma_gene %>%
  rename(magma_p = P) %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_logP = -log10(magma_p),
    magma_rank = rank(magma_p, ties.method = "min")
  ) %>%
  arrange(magma_p)

cat("MAGMA genes:", nrow(magma_gene), "\n")
cat("MAGMA genes with FDR < 0.05:", sum(magma_gene$magma_fdr < 0.05), "\n")

# Load Pathway results
pathway_results <- read.csv(
  file.path(OUT_DIR, "sorghum_stem_vol_catfish_results_MVN.csv"),
  stringsAsFactors = FALSE
)

# Bonferroni correction for pathways
n_pathways <- nrow(pathway_results)
pathway_results <- pathway_results %>%
  mutate(
    pathway_bonf = omni_p_final * n_pathways,
    pathway_bonf = pmin(pathway_bonf, 1),
    sig_bonf = pathway_bonf < 0.05
  ) %>%
  arrange(omni_p_final) %>%
  mutate(pathway_rank = row_number())

cat("Total pathways:", n_pathways, "\n")
cat("Pathways with Bonferroni < 0.05:", sum(pathway_results$sig_bonf), "\n")

# Get genes in significant pathways (Bonferroni)
sig_pathways <- pathway_results %>% filter(sig_bonf)
pathway_genes_list <- strsplit(sig_pathways$genes_used, ";")
names(pathway_genes_list) <- sig_pathways$pathway_id

pathway_genes <- data.frame(
  pathway_id = rep(sig_pathways$pathway_id, lengths(pathway_genes_list)),
  pathway_rank = rep(sig_pathways$pathway_rank, lengths(pathway_genes_list)),
  pathway_p = rep(sig_pathways$omni_p_final, lengths(pathway_genes_list)),
  GENE = unlist(pathway_genes_list),
  stringsAsFactors = FALSE
) %>%
  mutate(GENE = trimws(GENE))

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(
    n_pathways = n_distinct(pathway_id),
    best_pathway_rank = min(pathway_rank),
    best_pathway_p = min(pathway_p),
    pathways = paste(pathway_id, collapse = "; "),
    .groups = "drop"
  )

cat("Genes in Bonferroni-significant pathways:", nrow(pathway_gene_support), "\n")

# ==============================================================================
# Integrate all layers
# ==============================================================================

gene_evidence <- magma_gene %>%
  select(GENE, CHR, START, STOP, magma_p, magma_fdr, magma_logP, magma_rank) %>%
  left_join(gwas_gene, by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas = !is.na(gwas_logP) & gwas_logP > 7,
    hit_magma = !is.na(magma_fdr) & magma_fdr < 0.05,
    hit_pathway = !is.na(n_pathways) & n_pathways >= 1,
    in_pathway = !is.na(n_pathways)
  )

cat("\nHits summary:\n")
cat("  GWAS -log10(p) > 7:", sum(gene_evidence$hit_gwas, na.rm = TRUE), "\n")
cat("  MAGMA FDR < 0.05:", sum(gene_evidence$hit_magma, na.rm = TRUE), "\n")
cat("  Pathway Bonferroni < 0.05:", sum(gene_evidence$hit_pathway, na.rm = TRUE), "\n")

# ==============================================================================
# Panel A: UpSet Plot
# ==============================================================================

cat("\n=== Creating Panel A: UpSet Plot ===\n")

upset_df <- gene_evidence %>%
  transmute(
    GENE,
    `GWAS (-log10(p) > 7)` = as.integer(hit_gwas),
    `MAGMA (FDR < 0.05)` = as.integer(hit_magma),
    `Pathway (Bonf. < 0.05)` = as.integer(hit_pathway)
  ) %>%
  as.data.frame()

rownames(upset_df) <- upset_df$GENE
upset_df$GENE <- NULL
upset_df <- upset_df[rowSums(upset_df) > 0, ]

cat("Genes with at least one hit:", nrow(upset_df), "\n")

png(file.path(OUT_DIR, "Figure5_A_UpSet.png"),
    width = 4000, height = 3200, res = 300)
upset(
  upset_df,
  sets = c("GWAS (-log10(p) > 7)", "MAGMA (FDR < 0.05)", "Pathway (Bonf. < 0.05)"),
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(3, 3, 3, 3, 3, 3),
  point.size = 8,
  line.size = 2.5,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = c("coral", "steelblue", "forestgreen"),
  scale.sets = "identity",
  set_size.angles = 45
)
dev.off()

cat("Panel A saved\n")

# ==============================================================================
# Panel B: Summary Bar Chart
# ==============================================================================

cat("\n=== Creating Panel B: Summary Bar Chart ===\n")

summary_counts <- data.frame(
  Category = c(
    "All 3",
    "GWAS + Pathway",
    "MAGMA + Pathway",
    "GWAS + MAGMA",
    "GWAS only",
    "MAGMA only",
    "Pathway only"
  ),
  Count = c(
    sum(gene_evidence$hit_gwas & gene_evidence$hit_magma & gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & gene_evidence$hit_pathway & !gene_evidence$hit_magma, na.rm = TRUE),
    sum(gene_evidence$hit_magma & gene_evidence$hit_pathway & !gene_evidence$hit_gwas, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & gene_evidence$hit_magma & !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & !gene_evidence$hit_magma & !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_magma & !gene_evidence$hit_gwas & !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_pathway & !gene_evidence$hit_gwas & !gene_evidence$hit_magma, na.rm = TRUE)
  )
)

summary_counts$Category <- factor(summary_counts$Category, levels = summary_counts$Category)

panel_b <- ggplot(summary_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 8, fontface = "bold") +
  scale_fill_manual(values = c(
    "All 3" = "red",
    "GWAS + Pathway" = "hotpink",
    "MAGMA + Pathway" = "darkcyan",
    "GWAS + MAGMA" = "purple",
    "GWAS only" = "orange",
    "MAGMA only" = "navy",
    "Pathway only" = "forestgreen"
  )) +
  labs(
    x = "", y = "Number of Genes"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1, size = 20),
    legend.position = "none"
  ) +
  ylim(0, max(summary_counts$Count) * 1.15)

ggsave(file.path(OUT_DIR, "Figure5_B_Summary.png"), panel_b,
       width = 12, height = 10, dpi = 300, bg = "white")

cat("Panel B saved\n")

# ==============================================================================
# Panel C: Density Plot for GWAS and MAGMA
# ==============================================================================

cat("\n=== Creating Panel C: Density Plot ===\n")

# Prepare data for density plot
density_data <- gene_evidence %>%
  filter(!is.na(gwas_logP) & !is.na(magma_logP)) %>%
  mutate(Group = ifelse(in_pathway, "In Significant Pathways", "Not in Significant Pathways")) %>%
  select(GENE, Group, gwas_logP, magma_logP) %>%
  pivot_longer(cols = c(gwas_logP, magma_logP),
               names_to = "Method",
               values_to = "logP") %>%
  mutate(Method = ifelse(Method == "gwas_logP", "GWAS", "MAGMA"))

# Wilcoxon tests
wilcox_gwas <- wilcox.test(
  gene_evidence$gwas_logP[gene_evidence$in_pathway],
  gene_evidence$gwas_logP[!gene_evidence$in_pathway],
  alternative = "greater"
)

wilcox_magma <- wilcox.test(
  gene_evidence$magma_logP[gene_evidence$in_pathway],
  gene_evidence$magma_logP[!gene_evidence$in_pathway],
  alternative = "greater"
)

cat("Wilcoxon test (pathway genes have higher -log10(p)):\n")
cat("  GWAS p-value:", format(wilcox_gwas$p.value, digits = 4), "\n")
cat("  MAGMA p-value:", format(wilcox_magma$p.value, digits = 4), "\n")

panel_c <- ggplot(density_data, aes(x = logP, fill = Group)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Method, scales = "free_x", ncol = 2) +
  scale_fill_manual(values = c("In Significant Pathways" = "forestgreen",
                                "Not in Significant Pathways" = "gray50")) +
  labs(
    x = expression(-log[10](p)),
    y = "Density",
    caption = paste0("Wilcoxon p: GWAS = ", format(wilcox_gwas$p.value, digits = 3),
                     ", MAGMA = ", format(wilcox_magma$p.value, digits = 3))
  ) +
  plot_theme +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 20),
    strip.text = element_text(size = 24, face = "bold"),
    plot.caption = element_text(size = 18, hjust = 0.5)
  )

ggsave(file.path(OUT_DIR, "Figure5_C_Density.png"), panel_c,
       width = 14, height = 8, dpi = 300, bg = "white")

cat("Panel C saved\n")

# ==============================================================================
# Panel D: GWAS vs MAGMA Scatter
# ==============================================================================

cat("\n=== Creating Panel D: Scatter Plot ===\n")

scatter_df <- gene_evidence %>%
  filter(!is.na(gwas_logP) & !is.na(magma_logP)) %>%
  mutate(
    highlight = case_when(
      hit_gwas & hit_magma & hit_pathway ~ "All 3",
      hit_gwas & hit_pathway & !hit_magma ~ "GWAS + Pathway",
      hit_magma & hit_pathway & !hit_gwas ~ "MAGMA + Pathway",
      hit_gwas & hit_magma & !hit_pathway ~ "GWAS + MAGMA",
      hit_gwas & !hit_magma & !hit_pathway ~ "GWAS only",
      hit_magma & !hit_gwas & !hit_pathway ~ "MAGMA only",
      hit_pathway & !hit_gwas & !hit_magma ~ "Pathway only",
      TRUE ~ "None"
    ),
    highlight = factor(highlight, levels = c(
      "All 3", "GWAS + Pathway", "MAGMA + Pathway", "GWAS + MAGMA",
      "GWAS only", "MAGMA only", "Pathway only", "None"
    )),
    point_size = case_when(
      highlight == "All 3" ~ 6,
      highlight == "None" ~ 1.5,
      TRUE ~ 3.5
    )
  )

# Thresholds
gwas_threshold <- 7
magma_fdr_genes <- gene_evidence %>% filter(magma_fdr < 0.05)
magma_threshold <- min(magma_fdr_genes$magma_logP, na.rm = TRUE)

panel_d <- ggplot(scatter_df, aes(x = gwas_logP, y = magma_logP, color = highlight, size = point_size)) +
  geom_point(alpha = 0.6) +
  scale_size_identity() +
  geom_hline(yintercept = magma_threshold, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = gwas_threshold, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c(
    "All 3" = "red",
    "GWAS + Pathway" = "hotpink",
    "MAGMA + Pathway" = "darkcyan",
    "GWAS + MAGMA" = "purple",
    "GWAS only" = "orange",
    "MAGMA only" = "navy",
    "Pathway only" = "forestgreen",
    "None" = "gray80"
  )) +
  labs(
    x = expression(-log[10](GWAS~p)),
    y = expression(-log[10](MAGMA~p)),
    color = ""
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    legend.justification = "center",
    legend.text = element_text(size = 18)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE,
                               override.aes = list(size = 6, alpha = 1)))

ggsave(file.path(OUT_DIR, "Figure5_D_Scatter.png"), panel_d,
       width = 12, height = 10, dpi = 300, bg = "white")

cat("Panel D saved\n")

# ==============================================================================
# Combine all panels into Figure 5
# ==============================================================================

cat("\n=== Combining panels ===\n")

# Use cowplot for better image handling
# Create top row: UpSet (A) + Summary bar (B)
panel_a_plot <- ggdraw() +
  draw_image(file.path(OUT_DIR, "Figure5_A_UpSet.png")) +
  draw_label("A", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

# Add simple label to panel_b
panel_b_plot <- ggdraw(panel_b) +
  draw_label("B", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

# Add simple label to panel_c
panel_c_plot <- ggdraw(panel_c) +
  draw_label("C", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

# Add simple label to panel_d
panel_d_plot <- ggdraw(panel_d) +
  draw_label("D", x = 0.02, y = 0.98, hjust = 0, vjust = 1,
             fontface = "bold", size = 32)

# Combine using cowplot - top row and bottom row
top_row <- plot_grid(panel_a_plot, panel_b_plot, ncol = 2,
                     rel_widths = c(1.3, 1), align = "h")
bottom_row <- plot_grid(panel_c_plot, panel_d_plot, ncol = 2,
                        rel_widths = c(1.2, 1), align = "h")

fig5 <- plot_grid(top_row, bottom_row, nrow = 2,
                  rel_heights = c(1.2, 1))

ggsave(file.path(OUT_DIR, "Figure5_Candidate_Genes.png"), fig5,
       width = 24, height = 26, dpi = 300, bg = "white")

ggsave(file.path(OUT_DIR, "Figure5_Candidate_Genes.pdf"), fig5,
       width = 24, height = 26, bg = "white")

cat("Figure 5 saved!\n")

# ==============================================================================
# Final Table: Genes hit in all three layers
# ==============================================================================

cat("\n=== Creating Final Table ===\n")

# Get genes that are hits in all three
all_three_genes <- gene_evidence %>%
  filter(hit_gwas & hit_magma & hit_pathway) %>%
  select(
    GENE, CHR, START, STOP,
    gwas_min_p, gwas_logP, gwas_rank,
    magma_p, magma_fdr, magma_logP, magma_rank,
    n_pathways, best_pathway_rank, best_pathway_p, pathways
  ) %>%
  arrange(magma_rank)

cat("Genes in all three layers:", nrow(all_three_genes), "\n")

# Save the table
write.csv(
  all_three_genes,
  file.path(OUT_DIR, "Table_Genes_All_Three_Layers.csv"),
  row.names = FALSE
)

# Print summary
cat("\n=== Genes in All Three Layers ===\n")
cat(strrep("-", 80), "\n")
for (i in seq_len(nrow(all_three_genes))) {
  g <- all_three_genes[i, ]
  cat(sprintf("%2d. %s (Chr%s: %d-%d)\n", i, g$GENE, g$CHR, g$START, g$STOP))
  cat(sprintf("    GWAS: -log10(p) = %.2f, rank = %d\n", g$gwas_logP, g$gwas_rank))
  cat(sprintf("    MAGMA: -log10(p) = %.2f, FDR = %.2e, rank = %d\n",
              g$magma_logP, g$magma_fdr, g$magma_rank))
  cat(sprintf("    Pathways (%d): rank = %d, %s\n",
              g$n_pathways, g$best_pathway_rank, g$pathways))
  cat("\n")
}

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== DONE ===\n")
cat("Figure 5 saved:\n")
cat("  -", file.path(OUT_DIR, "Figure5_Candidate_Genes.png"), "\n")
cat("  -", file.path(OUT_DIR, "Figure5_Candidate_Genes.pdf"), "\n")
cat("\nIndividual panels:\n")
cat("  -", file.path(OUT_DIR, "Figure5_A_UpSet.png"), "\n")
cat("  -", file.path(OUT_DIR, "Figure5_B_Summary.png"), "\n")
cat("  -", file.path(OUT_DIR, "Figure5_C_Density.png"), "\n")
cat("  -", file.path(OUT_DIR, "Figure5_D_Scatter.png"), "\n")
cat("\nTable saved:\n")
cat("  -", file.path(OUT_DIR, "Table_Genes_All_Three_Layers.csv"), "-", nrow(all_three_genes), "genes\n")
