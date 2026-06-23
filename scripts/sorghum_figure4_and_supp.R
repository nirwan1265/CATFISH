# ==============================================================================
# Figure 4: GWAS, MAGMA, and Pathway Manhattan Plots
# Sorghum Stem Volume
# ==============================================================================

library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)

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
# Plot Theme
# ==============================================================================

plot_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Chromosome colors (alternating)
chr_colors <- rep(c("gray30", "gray60"), 5)

# ==============================================================================
# Panel A: GWAS Manhattan Plot
# ==============================================================================

cat("=== Loading GWAS data ===\n")

gwas_raw <- fread(

  "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt",
  header = TRUE
)

names(gwas_raw)[names(gwas_raw) == "p_wald"] <- "P"
names(gwas_raw)[names(gwas_raw) == "ps"] <- "POS"
names(gwas_raw)[names(gwas_raw) == "chr"] <- "CHR"

cat("Total SNPs:", nrow(gwas_raw), "\n")

# Calculate cumulative positions for Manhattan plot
gwas_raw <- gwas_raw %>%
  filter(!is.na(P) & P > 0) %>%
  mutate(logP = -log10(P))

# Get chromosome lengths
chr_lengths <- gwas_raw %>%
  group_by(CHR) %>%
  summarise(chr_len = max(POS), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    chr_mid = cum_start + chr_len / 2
  )

gwas_plot_data <- gwas_raw %>%
  left_join(chr_lengths %>% select(CHR, cum_start), by = "CHR") %>%
  mutate(cum_pos = POS + cum_start)

# Downsample for plotting (keep all significant, sample non-significant)
set.seed(42)
gwas_sig <- gwas_plot_data %>% filter(logP >= 5)
gwas_nonsig <- gwas_plot_data %>% filter(logP < 5) %>% sample_frac(0.05)
gwas_plot_sample <- bind_rows(gwas_sig, gwas_nonsig)

cat("Plotting", nrow(gwas_plot_sample), "SNPs\n")

# Genome-wide significance threshold
gwas_threshold <- -log10(5e-8)
suggestive_threshold <- -log10(1e-5)

panel_a <- ggplot(gwas_plot_sample, aes(x = cum_pos, y = logP, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = gwas_threshold, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = suggestive_threshold, color = "blue", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60")) +
  scale_x_continuous(
    breaks = chr_lengths$chr_mid,
    labels = chr_lengths$CHR,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title = "A. GWAS Manhattan Plot",
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  plot_theme +
  coord_cartesian(ylim = c(0, max(gwas_plot_sample$logP) * 1.05))

cat("Panel A done\n")

# ==============================================================================
# Panel B: MAGMA Gene Manhattan Plot
# ==============================================================================

cat("\n=== Loading MAGMA data ===\n")

magma_gene <- read.table(
  file.path(OUT_DIR, "sorghum_stem_vol_genes_combined.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

cat("Total genes:", nrow(magma_gene), "\n")

magma_gene <- magma_gene %>%
  mutate(
    logP = -log10(P),
    gene_mid = (START + STOP) / 2
  )

# Calculate cumulative positions
chr_lengths_magma <- magma_gene %>%
  group_by(CHR) %>%
  summarise(chr_len = max(STOP), .groups = "drop") %>%
  arrange(CHR) %>%
  mutate(
    cum_start = lag(cumsum(chr_len), default = 0),
    chr_mid = cum_start + chr_len / 2
  )

magma_plot_data <- magma_gene %>%
  left_join(chr_lengths_magma %>% select(CHR, cum_start), by = "CHR") %>%
  mutate(cum_pos = gene_mid + cum_start)

# MAGMA threshold (Bonferroni)
magma_bonf_threshold <- -log10(0.05 / nrow(magma_gene))
magma_fdr_genes <- magma_gene %>%
  mutate(fdr = p.adjust(P, method = "BH")) %>%
  filter(fdr < 0.05)
magma_fdr_threshold <- -log10(max(magma_fdr_genes$P))

panel_b <- ggplot(magma_plot_data, aes(x = cum_pos, y = logP, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = magma_bonf_threshold, color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = magma_fdr_threshold, color = "blue", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60")) +
  scale_x_continuous(
    breaks = chr_lengths_magma$chr_mid,
    labels = chr_lengths_magma$CHR,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title = "B. MAGMA Gene Manhattan Plot",
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  plot_theme +
  coord_cartesian(ylim = c(0, max(magma_plot_data$logP) * 1.05))

cat("Panel B done\n")

# ==============================================================================
# Panel C: CATFISH Pathway Results - Component Tests Heatmap
# ==============================================================================

cat("\n=== Loading Pathway data ===\n")
cat("Using:", PATHWAY_FILE, "\n")

pathway_results <- read.csv(PATHWAY_FILE, stringsAsFactors = FALSE)

cat("Total pathways:", nrow(pathway_results), "\n")

# Select top 20 pathways for detailed view
top_pathways <- pathway_results %>%
  arrange(omni_p_final) %>%
  head(20) %>%
  mutate(
    pathway_label = paste0(pathway_name, " (", n_genes, ")"),
    pathway_label = factor(pathway_label, levels = rev(pathway_label))
  )

# Prepare component test data for heatmap
component_data <- top_pathways %>%
  select(pathway_label, n_genes,
         ACAT = acat_p,
         Fisher = fisher_p,
         Stouffer = stouffer_p_analytic,
         minP = minp_p_analytic,
         TFisher = tfisher_p_analytic,
         Omnibus = omni_p_final) %>%
  tidyr::pivot_longer(
    cols = c(ACAT, Fisher, Stouffer, minP, TFisher, Omnibus),
    names_to = "Method",
    values_to = "p_value"
  ) %>%
  mutate(
    logP = -log10(pmax(p_value, 1e-16)),
    Method = factor(Method, levels = c("ACAT", "Fisher", "Stouffer", "minP", "TFisher", "Omnibus"))
  )

# Panel C: Component test heatmap only
panel_c <- ggplot(component_data, aes(x = Method, y = pathway_label, fill = logP)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "gray95", mid = "steelblue", high = "darkred",
    midpoint = 5, limits = c(0, 16),
    name = expression(-log[10](p))
  ) +
  labs(
    title = "C. CATFISH Pathway Analysis",
    x = "",
    y = ""
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, vjust = 1, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 10),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

cat("Panel C done\n")

# ==============================================================================
# Combine Panels
# ==============================================================================

cat("\n=== Combining panels ===\n")

# Combine A, B, C in one figure
fig4 <- (panel_a / panel_b / panel_c) +
  plot_layout(heights = c(1, 1, 1.2))

ggsave(
  file.path(OUT_DIR, paste0("Figure4_GWAS_MAGMA_Pathway", OUT_SUFFIX, ".png")),
  fig4,
  width = 12, height = 16,
  dpi = 300, bg = "white"
)

ggsave(
  file.path(OUT_DIR, paste0("Figure4_GWAS_MAGMA_Pathway", OUT_SUFFIX, ".pdf")),
  fig4,
  width = 12, height = 16,
  bg = "white"
)

cat("Figure 4 saved!\n")

# ==============================================================================
# Supplementary Table: Significant Results
# ==============================================================================

cat("\n=== Creating Supplementary Tables ===\n")

# --- Table S1: GWAS genes with -log10(p) > 7 ---
gene_loc <- read.delim("inst/extdata/sorghum.genes.loc", stringsAsFactors = FALSE)
gene_loc <- gene_loc %>%
  mutate(
    START_EXT = pmax(0, START - 25000),
    STOP_EXT = STOP + 25000
  )

# Map GWAS SNPs to genes
gwas_dt <- as.data.table(gwas_raw)
gene_dt <- as.data.table(gene_loc)
setkey(gene_dt, CHR, START_EXT, STOP_EXT)
gwas_dt[, c("START_EXT", "STOP_EXT") := .(POS, POS)]
setkey(gwas_dt, CHR, START_EXT, STOP_EXT)

snp_gene_map <- foverlaps(
  gwas_dt,
  gene_dt,
  by.x = c("CHR", "START_EXT", "STOP_EXT"),
  by.y = c("CHR", "START_EXT", "STOP_EXT"),
  type = "within",
  nomatch = NULL
)

snp_gene_df <- as.data.frame(snp_gene_map)
# Fix column names after foverlaps
if ("i.CHR" %in% names(snp_gene_df)) {
  names(snp_gene_df)[names(snp_gene_df) == "i.CHR"] <- "SNP_CHR"
}

gwas_gene_summary <- snp_gene_df %>%
  group_by(GENE) %>%
  summarise(
    CHR = first(CHR),
    gene_start = first(START),
    gene_stop = first(STOP),
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_logP = -log10(min(P, na.rm = TRUE)),
    n_snps_mapped = n(),
    .groups = "drop"
  ) %>%
  filter(gwas_logP > 7) %>%
  arrange(desc(gwas_logP))

cat("GWAS genes with -log10(p) > 7:", nrow(gwas_gene_summary), "\n")

write.csv(
  gwas_gene_summary,
  file.path(OUT_DIR, "TableS1_GWAS_genes_logP_gt7.csv"),
  row.names = FALSE
)

# --- Table S2: MAGMA genes with FDR < 0.05 ---
magma_sig <- magma_gene %>%
  mutate(
    FDR = p.adjust(P, method = "BH"),
    Bonferroni = p.adjust(P, method = "bonferroni")
  ) %>%
  filter(FDR < 0.05) %>%
  select(GENE, CHR, START, STOP, NSNPS, ZSTAT, P, FDR, Bonferroni) %>%
  arrange(P)

cat("MAGMA genes with FDR < 0.05:", nrow(magma_sig), "\n")

write.csv(
  magma_sig,
  file.path(OUT_DIR, "TableS2_MAGMA_genes_FDR_lt0.05.csv"),
  row.names = FALSE
)

# --- Table S3: Pathways with BH < 0.05 and their genes ---
pathway_sig <- pathway_results %>%
  filter(omni_p_final_BH < 0.05) %>%
  select(pathway_id, pathway_name, n_genes,
         acat_p, fisher_p, stouffer_p_analytic, minp_p_analytic, tfisher_p_analytic,
         omni_p_final, FDR = omni_p_final_BH, genes_used) %>%
  arrange(omni_p_final)

cat("Pathways with FDR < 0.05:", nrow(pathway_sig), "\n")

write.csv(
  pathway_sig,
  file.path(OUT_DIR, paste0("TableS3_Pathways_FDR_lt0.05", OUT_SUFFIX, ".csv")),
  row.names = FALSE
)

# --- Table S4: All genes in significant pathways ---
pathway_genes_expanded <- pathway_sig %>%
  select(pathway_id, pathway_name, omni_p_final, pathway_FDR = FDR, genes_used) %>%
  mutate(genes_list = strsplit(genes_used, ";")) %>%
  tidyr::unnest(genes_list) %>%
  mutate(gene = trimws(genes_list)) %>%
  select(pathway_id, pathway_name, pathway_omni_p = omni_p_final, pathway_FDR, gene)

# Add MAGMA and GWAS info for each gene
magma_for_join <- magma_gene %>%
  mutate(magma_FDR = p.adjust(P, method = "BH")) %>%
  select(GENE, magma_p = P, magma_FDR)

gwas_for_join <- gwas_gene_summary %>%
  select(GENE, gwas_min_p, gwas_logP)

pathway_genes_full <- pathway_genes_expanded %>%
  left_join(magma_for_join, by = c("gene" = "GENE")) %>%
  left_join(gwas_for_join, by = c("gene" = "GENE")) %>%
  arrange(pathway_FDR, pathway_id, magma_p)

cat("Total gene-pathway associations:", nrow(pathway_genes_full), "\n")

write.csv(
  pathway_genes_full,
  file.path(OUT_DIR, paste0("TableS4_Pathway_genes_detailed", OUT_SUFFIX, ".csv")),
  row.names = FALSE
)

# ==============================================================================
# Summary
# ==============================================================================

cat("\n=== DONE ===\n")
cat("Figure saved:\n")
cat("  -", file.path(OUT_DIR, "Figure4_GWAS_MAGMA_Pathway.png"), "\n")
cat("  -", file.path(OUT_DIR, "Figure4_GWAS_MAGMA_Pathway.pdf"), "\n")
cat("\nSupplementary tables saved:\n")
cat("  -", file.path(OUT_DIR, "TableS1_GWAS_genes_logP_gt7.csv"), "-", nrow(gwas_gene_summary), "genes\n")
cat("  -", file.path(OUT_DIR, "TableS2_MAGMA_genes_FDR_lt0.05.csv"), "-", nrow(magma_sig), "genes\n")
cat("  -", file.path(OUT_DIR, "TableS3_Pathways_FDR_lt0.05.csv"), "-", nrow(pathway_sig), "pathways\n")
cat("  -", file.path(OUT_DIR, "TableS4_Pathway_genes_detailed.csv"), "-", nrow(pathway_genes_full), "gene-pathway pairs\n")
