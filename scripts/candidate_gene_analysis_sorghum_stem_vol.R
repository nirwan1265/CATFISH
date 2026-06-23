# ==============================================================================
# Candidate Gene Analysis: Integrating GWAS, MAGMA, and Pathway Evidence
# Sorghum - Stem Volume Phenotype
# ==============================================================================
#
# Three layers of evidence:
#   1. GWAS: SNP-level association (mapped to genes via 25kb window)
#   2. MAGMA: Gene-level association (LD-aware multi-SNP aggregation)
#   3. Pathway: Set-level association (coordinated/diffuse signal detection)
#
# ==============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(data.table)

if (!requireNamespace("UpSetR", quietly = TRUE)) {
  install.packages("UpSetR")
}
library(UpSetR)

if (!requireNamespace("patchwork", quietly = TRUE)) {
  install.packages("patchwork")
}
library(patchwork)

# ==============================================================================
# USER PARAMETERS
# ==============================================================================

# ---------- GWAS layer ----------
GWAS_MODE           <- "top_pct"   # "threshold" or "top_pct"
GWAS_P_THRESHOLD    <- 1e-5        # p-value cutoff (GWAS_MODE="threshold")
GWAS_TOP_PCT        <- 1           # top X% of genes (GWAS_MODE="top_pct")

# ---------- MAGMA layer ----------
MAGMA_MODE          <- "top_pct"   # "threshold" or "top_pct"
MAGMA_FDR_THRESHOLD <- 0.3         # FDR cutoff (MAGMA_MODE="threshold")
MAGMA_TOP_PCT       <- 5           # top X% of genes (MAGMA_MODE="top_pct")

# ---------- Pathway layer ----------
PATHWAY_MODE          <- "fdr"     # "top_k" or "fdr"
TOP_K_PATHWAYS        <- 20        # count   (PATHWAY_MODE="top_k")
PATHWAY_FDR_THRESHOLD <- 0.05      # FDR cutoff (PATHWAY_MODE="fdr")

# ---------- General ----------
WINDOW_SIZE         <- 25000       # SNP-to-gene mapping window (bp)
OUT_DIR             <- "sorghum_catfish_results"

# ==============================================================================
# 0. Plot Theme
# ==============================================================================

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 24,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(size = 24, face = "bold"),
    axis.title.y   = element_text(size = 24, face = "bold"),
    axis.text.x    = element_text(size = 24, color = "black"),
    axis.text.y    = element_text(size = 24, color = "black"),
    axis.line      = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 24)
  )

# ==============================================================================
# 1. Load Gene Location File
# ==============================================================================

gene_loc <- read.delim(
  "inst/extdata/sorghum.genes.loc",
  stringsAsFactors = FALSE
)

head(gene_loc)
cat("Total genes in annotation:", nrow(gene_loc), "\n")

gene_loc <- gene_loc %>%
  mutate(
    START_EXT = pmax(0, START - WINDOW_SIZE),
    STOP_EXT = STOP + WINDOW_SIZE
  )

# ==============================================================================
# 2. Load GWAS Results and Map SNPs to Genes
# ==============================================================================

gwas_raw <- fread(
  "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt",
  header = TRUE
)

cat("Total SNPs in GWAS:", nrow(gwas_raw), "\n")

# Standardize column names (GEMMA format: chr, rs, ps, p_wald)
if ("p_wald" %in% names(gwas_raw)) {
  names(gwas_raw)[names(gwas_raw) == "p_wald"] <- "P"
}
if ("ps" %in% names(gwas_raw)) {
  names(gwas_raw)[names(gwas_raw) == "ps"] <- "POS"
}
if ("chr" %in% names(gwas_raw)) {
  names(gwas_raw)[names(gwas_raw) == "chr"] <- "CHR"
}
if ("rs" %in% names(gwas_raw)) {
  names(gwas_raw)[names(gwas_raw) == "rs"] <- "SNP_ID"
}

cat("\nGWAS p-value summary:\n")
print(summary(gwas_raw$P))

cat("\nSNPs at different thresholds:\n")
cat("  p < 1e-7:", sum(gwas_raw$P < 1e-7, na.rm = TRUE), "\n")
cat("  p < 1e-6:", sum(gwas_raw$P < 1e-6, na.rm = TRUE), "\n")
cat("  p < 1e-5:", sum(gwas_raw$P < 1e-5, na.rm = TRUE), "\n")
cat("  p < 1e-4:", sum(gwas_raw$P < 1e-4, na.rm = TRUE), "\n")
cat("  p < 1e-3:", sum(gwas_raw$P < 1e-3, na.rm = TRUE), "\n")
cat("  p < 0.05:", sum(gwas_raw$P < 0.05, na.rm = TRUE), "\n")

cat("\nMapping SNPs to genes (25kb window)...\n")

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

cat("SNP-gene mappings:", nrow(snp_gene_map), "\n")
cat("Unique genes with mapped SNPs:", length(unique(snp_gene_map$GENE)), "\n")

gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_n_snps = n(),
    .groups = "drop"
  ) %>%
  arrange(gwas_min_p) %>%
  filter(!duplicated(GENE)) %>%
  mutate(gwas_rank = rank(gwas_min_p, ties.method = "min"))

cat("\nGenes with GWAS evidence:", nrow(gwas_gene), "\n")

cat("\nGenes with GWAS hit (min SNP p-value):\n")
cat("  p < 1e-7:", sum(gwas_gene$gwas_min_p < 1e-7), "\n")
cat("  p < 1e-6:", sum(gwas_gene$gwas_min_p < 1e-6), "\n")
cat("  p < 1e-5:", sum(gwas_gene$gwas_min_p < 1e-5), "\n")
cat("  p < 1e-4:", sum(gwas_gene$gwas_min_p < 1e-4), "\n")
cat("  p < 1e-3:", sum(gwas_gene$gwas_min_p < 1e-3), "\n")

# Compute effective GWAS threshold based on mode
if (GWAS_MODE == "top_pct") {
  gwas_n_top <- ceiling(nrow(gwas_gene) * GWAS_TOP_PCT / 100)
  gwas_eff_threshold <- sort(gwas_gene$gwas_min_p)[gwas_n_top]
  cat("\nGWAS mode: top", GWAS_TOP_PCT, "%\n")
  cat("  Effective p cutoff:",
      format(gwas_eff_threshold, digits = 3), "\n")
  cat("  Genes selected:", gwas_n_top,
      "of", nrow(gwas_gene), "\n")
} else {
  gwas_eff_threshold <- GWAS_P_THRESHOLD
  cat("\nGWAS mode: p-value threshold <",
      GWAS_P_THRESHOLD, "\n")
  cat("  Genes selected:",
      sum(gwas_gene$gwas_min_p < gwas_eff_threshold), "\n")
}

# ==============================================================================
# 3. Load MAGMA Gene Results
# ==============================================================================

magma_gene <- read.table(
  file.path(OUT_DIR, "sorghum_stem_vol_genes_combined.txt"),
  header = TRUE,
  stringsAsFactors = FALSE
)

cat("\nMAGMA columns:", paste(names(magma_gene), collapse = ", "), "\n")
cat("Genes with MAGMA results:", nrow(magma_gene), "\n")

# Use the P column (already created during MAGMA processing)
magma_gene <- magma_gene %>%
  rename(magma_p = P) %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_bonf = p.adjust(magma_p, method = "bonferroni")
  ) %>%
  arrange(magma_p)

cat("\nMAGMA significant genes:\n")
cat("  Bonferroni < 0.05:", sum(magma_gene$magma_bonf < 0.05, na.rm = TRUE), "\n")
cat("  FDR < 0.05:", sum(magma_gene$magma_fdr < 0.05, na.rm = TRUE), "\n")
cat("  FDR < 0.10:", sum(magma_gene$magma_fdr < 0.10, na.rm = TRUE), "\n")
cat("  Nominal p < 0.05:", sum(magma_gene$magma_p < 0.05, na.rm = TRUE), "\n")

cat("\nTop 10 MAGMA genes:\n")
print(head(magma_gene, 10))

# Compute effective MAGMA threshold based on mode
if (MAGMA_MODE == "top_pct") {
  magma_n_top <- ceiling(nrow(magma_gene) * MAGMA_TOP_PCT / 100)
  magma_eff_threshold <- sort(magma_gene$magma_p)[magma_n_top]
  cat("\nMAGMA mode: top", MAGMA_TOP_PCT, "%\n")
  cat("  Effective p cutoff:",
      format(magma_eff_threshold, digits = 3), "\n")
  cat("  Genes selected:", magma_n_top,
      "of", nrow(magma_gene), "\n")
} else {
  magma_eff_threshold <- NULL
  cat("\nMAGMA mode: FDR threshold <",
      MAGMA_FDR_THRESHOLD, "\n")
  cat("  Genes selected:",
      sum(magma_gene$magma_fdr < MAGMA_FDR_THRESHOLD),
      "\n")
}

# ==============================================================================
# 4. Load Pathway Results
# ==============================================================================

omni_results <- read.csv(
  file.path(OUT_DIR, "sorghum_stem_vol_catfish_results_MVN.csv"),
  stringsAsFactors = FALSE
)

cat("\nTotal pathways:", nrow(omni_results), "\n")

omni_col <- "omni_p_final"
if (!omni_col %in% names(omni_results)) {
  omni_col <- "omni_p_mvn"
  if (!omni_col %in% names(omni_results)) {
    omni_col <- "omni_p_analytic"
  }
}
cat("Using pathway p-value column:", omni_col, "\n")

cat("\nPathway significant counts:\n")
bonf_thresh <- 0.05 / nrow(omni_results)
cat("  Bonferroni (", format(bonf_thresh, digits = 3), "):",
    sum(omni_results[[omni_col]] < bonf_thresh, na.rm = TRUE), "\n")
cat("  FDR < 0.10:",
    sum(p.adjust(omni_results[[omni_col]], "BH") < 0.10, na.rm = TRUE), "\n")
cat("  Nominal < 0.05:",
    sum(omni_results[[omni_col]] < 0.05, na.rm = TRUE), "\n")

# Select pathways based on mode
if (PATHWAY_MODE == "fdr") {
  omni_results$pathway_fdr <- p.adjust(
    omni_results[[omni_col]], "BH"
  )
  top_pathways <- omni_results %>%
    filter(pathway_fdr < PATHWAY_FDR_THRESHOLD) %>%
    arrange(!!sym(omni_col))
  cat("\nPathway mode: FDR <", PATHWAY_FDR_THRESHOLD, "\n")
  cat("  Significant pathways:", nrow(top_pathways), "\n")
} else {
  top_pathways <- omni_results %>%
    arrange(!!sym(omni_col)) %>%
    slice_head(n = TOP_K_PATHWAYS)
  cat("\nPathway mode: top", TOP_K_PATHWAYS, "\n")
}

n_pathways_used <- nrow(top_pathways)
cat("\nSelected", n_pathways_used, "pathways:\n")
for (i in seq_len(min(20, nrow(top_pathways)))) {
  cat("  ", i, ". ", top_pathways$pathway_id[i], " (",
      top_pathways$pathway_name[i], ")\n",
      "     p = ", format(top_pathways[[omni_col]][i], digits = 3), "\n",
      sep = "")
}

pathway_genes_list <- strsplit(top_pathways$genes_used, ";")
names(pathway_genes_list) <- top_pathways$pathway_id

pathway_genes <- data.frame(
  pathway_id = rep(top_pathways$pathway_id, lengths(pathway_genes_list)),
  pathway_name = rep(top_pathways$pathway_name, lengths(pathway_genes_list)),
  pathway_p = rep(top_pathways[[omni_col]], lengths(pathway_genes_list)),
  GENE = unlist(pathway_genes_list),
  stringsAsFactors = FALSE
)

pathway_genes$GENE <- trimws(pathway_genes$GENE)

cat("\nGenes in top pathways:", nrow(pathway_genes), "\n")
cat("Unique genes in top pathways:", length(unique(pathway_genes$GENE)), "\n")

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(
    n_top_pathways = n_distinct(pathway_id),
    best_pathway_p = min(pathway_p, na.rm = TRUE),
    mean_pathway_mlog10p = mean(-log10(pathway_p), na.rm = TRUE),
    pathways = paste(pathway_id, collapse = "; "),
    .groups = "drop"
  )

cat("Genes with pathway support:", nrow(pathway_gene_support), "\n")

# ==============================================================================
# 5. Integrate All Three Layers
# ==============================================================================

cat("\n\n========================================\n")
cat("Integrating Evidence Across Layers (Sorghum Stem Volume)\n")
cat("========================================\n")

gene_universe <- unique(magma_gene$GENE)
cat("\nGene universe (MAGMA-tested genes):", length(gene_universe), "\n")

gene_evidence <- magma_gene %>%
  left_join(gwas_gene, by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas = !is.na(gwas_min_p) &
      gwas_min_p <= gwas_eff_threshold,
    hit_magma = if (MAGMA_MODE == "top_pct") {
      !is.na(magma_p) & magma_p <= magma_eff_threshold
    } else {
      !is.na(magma_fdr) & magma_fdr < MAGMA_FDR_THRESHOLD
    },
    hit_pathway = !is.na(n_top_pathways) & n_top_pathways >= 1,
    support_layers = hit_gwas + hit_magma + hit_pathway,
    score = (hit_gwas * 1) + (hit_magma * 1) + (hit_pathway * 1) +
      0.2 * ifelse(!is.na(magma_p), -log10(magma_p), 0) +
      0.1 * ifelse(!is.na(gwas_min_p), -log10(gwas_min_p), 0) +
      0.1 * ifelse(!is.na(best_pathway_p), -log10(best_pathway_p), 0),
    gwas_rank = ifelse(is.na(gwas_min_p), NA_integer_, min_rank(gwas_min_p)),
    magma_rank = ifelse(is.na(magma_p), NA_integer_, min_rank(magma_p)),
    pathway_rank = ifelse(is.na(best_pathway_p), NA_integer_,
                          min_rank(best_pathway_p))
  ) %>%
  arrange(desc(score))

# ==============================================================================
# 6. UpSet Plot
# ==============================================================================

gwas_label <- if (GWAS_MODE == "top_pct") {
  paste0("GWAS (top ", GWAS_TOP_PCT, "%)")
} else {
  paste0("GWAS (p<", GWAS_P_THRESHOLD, ")")
}
magma_label <- if (MAGMA_MODE == "top_pct") {
  paste0("MAGMA (top ", MAGMA_TOP_PCT, "%)")
} else {
  paste0("MAGMA (FDR<", MAGMA_FDR_THRESHOLD, ")")
}
pathway_label <- if (PATHWAY_MODE == "fdr") {
  paste0("Pathway (FDR<", PATHWAY_FDR_THRESHOLD, ")")
} else {
  paste0("Pathway (top ", TOP_K_PATHWAYS, ")")
}

upset_df <- gene_evidence %>%
  transmute(
    GENE,
    !!gwas_label    := as.integer(hit_gwas),
    !!magma_label   := as.integer(hit_magma),
    !!pathway_label := as.integer(hit_pathway)
  ) %>%
  as.data.frame()

rownames(upset_df) <- upset_df$GENE
upset_df$GENE <- NULL

upset_df <- upset_df[rowSums(upset_df) > 0, ]

cat("\nGenes with at least one layer of evidence:", nrow(upset_df), "\n")

png(file.path(OUT_DIR, "candidate_genes_upset_sorghum_stem_vol.png"),
    width = 2100, height = 2000, res = 300)
upset(
  upset_df,
  sets = c(gwas_label, magma_label, pathway_label),
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = 2.2,
  point.size = 4.5,
  line.size = 1.5,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = c("coral", "steelblue", "forestgreen")
)
dev.off()

cat("UpSet plot saved to", file.path(OUT_DIR, "candidate_genes_upset_sorghum_stem_vol.png"), "\n")

# ==============================================================================
# 7. Top Candidate Genes Table
# ==============================================================================

cat("\n\n========================================\n")
cat("TOP CANDIDATE GENES (Sorghum Stem Volume)\n")
cat("========================================\n\n")

cat("Top 20 genes by composite score:\n")
cat(strrep("-", 80), "\n")

top20 <- head(gene_evidence, 20)

for (i in seq_len(nrow(top20))) {
  g <- top20[i, ]
  cat(
    sprintf("%2d. %s\n", i, g$GENE),
    sprintf("    Score: %.2f | Layers: %d\n", g$score, g$support_layers),
    sprintf("    GWAS min p: %s | MAGMA p: %s (FDR: %s)\n",
            ifelse(is.na(g$gwas_min_p), "NA",
                   format(g$gwas_min_p, digits = 3, scientific = TRUE)),
            format(g$magma_p, digits = 3, scientific = TRUE),
            format(g$magma_fdr, digits = 3)),
    sprintf("    Pathway: %s (n=%s)\n",
            ifelse(is.na(g$pathways), "None", g$pathways),
            ifelse(is.na(g$n_top_pathways), 0, g$n_top_pathways)),
    "\n", sep = ""
  )
}

write.csv(
  gene_evidence %>%
    select(GENE, magma_p, magma_fdr, magma_rank,
           gwas_min_p, gwas_n_snps, gwas_rank,
           n_top_pathways, best_pathway_p, pathway_rank, pathways,
           hit_gwas, hit_magma, hit_pathway,
           support_layers, score) %>%
    head(200),
  file.path(OUT_DIR, "candidate_genes_top200_sorghum_stem_vol.csv"),
  row.names = FALSE
)

cat("\nFull table saved to", file.path(OUT_DIR, "candidate_genes_top200_sorghum_stem_vol.csv"), "\n")

# ==============================================================================
# 8. Pathway Gene Enrichment Analysis
# ==============================================================================

cat("\n\n========================================\n")
cat("Pathway Gene Enrichment Analysis\n")
cat("========================================\n")

gene_evidence$in_pathway <- !is.na(gene_evidence$n_top_pathways)

pathway_genes_magma <- gene_evidence$magma_p[gene_evidence$in_pathway]
nonpathway_genes_magma <- gene_evidence$magma_p[!gene_evidence$in_pathway]

cat("Genes in top pathways:", length(pathway_genes_magma), "\n")
cat("Genes not in top pathways:", length(nonpathway_genes_magma), "\n")

wilcox_result <- wilcox.test(
  pathway_genes_magma,
  nonpathway_genes_magma,
  alternative = "less"
)

cat("\nWilcoxon rank-sum test (pathway genes have lower MAGMA p?):\n")
cat("  W =", wilcox_result$statistic, "\n")
cat("  p-value =", format(wilcox_result$p.value, digits = 4), "\n")

cat("\nMedian -log10(MAGMA p):\n")
cat("  Pathway genes:",
    round(median(-log10(pathway_genes_magma), na.rm = TRUE), 3), "\n")
cat("  Non-pathway genes:",
    round(median(-log10(nonpathway_genes_magma), na.rm = TRUE), 3), "\n")

density_df <- gene_evidence %>%
  mutate(Group = ifelse(in_pathway, "In Top Pathways", "Not in Top Pathways"))

density_plot <- ggplot(density_df, aes(x = -log10(magma_p), fill = Group)) +
  geom_density(alpha = 0.5) +
  labs(
    title = "MAGMA P-value Distribution (Sorghum Stem Volume)",
    subtitle = paste0("Wilcoxon p = ", format(wilcox_result$p.value, digits = 3)),
    x = expression(-log[10](MAGMA~p)),
    y = "Density"
  ) +
  scale_fill_manual(values = c("forestgreen", "gray50")) +
  theme_minimal() +
  theme(legend.position = "bottom") + plot_theme

ggsave(file.path(OUT_DIR, "pathway_gene_enrichment_density_sorghum_stem_vol.png"), density_plot,
       width = 8, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# 9. Scatter Plot: GWAS vs MAGMA
# ==============================================================================

# Calculate MAGMA p-value threshold for plotting (handle FDR mode)
if (MAGMA_MODE == "top_pct") {
  magma_p_threshold <- magma_eff_threshold
} else {
  # For FDR mode, find p-value corresponding to FDR threshold
  magma_p_threshold <- max(
    gene_evidence$magma_p[gene_evidence$magma_fdr < MAGMA_FDR_THRESHOLD],
    na.rm = TRUE
  )
}

scatter_df <- gene_evidence %>%
  filter(!is.na(gwas_min_p) & !is.na(magma_p)) %>%
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
    point_size = ifelse(highlight == "None", 2, 4)
  )

scatter_plot <- ggplot(
  scatter_df,
  aes(x = -log10(gwas_min_p), y = -log10(magma_p),
      color = highlight, size = point_size)
) +
  geom_point(alpha = 0.6) +
  scale_size_identity() +
  geom_hline(
    yintercept = -log10(magma_p_threshold),
    linetype = "dashed", color = "gray50"
  ) +
  geom_vline(
    xintercept = -log10(gwas_eff_threshold),
    linetype = "dashed", color = "gray50"
  ) +
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
    title = "GWAS vs MAGMA Evidence (Sorghum Stem Volume)",
    x = expression(-log[10](GWAS~min~p)),
    y = expression(-log[10](MAGMA~p)),
    color = "Evidence"
  ) +
  theme_minimal() +
  theme(legend.position = "right") + plot_theme +
  guides(
    color = guide_legend(
      nrow = 4, ncol = 2, byrow = TRUE,
      override.aes = list(size = 4, alpha = 1)
    )
  )

ggsave(file.path(OUT_DIR, "gwas_vs_magma_scatter_sorghum_stem_vol.png"), scatter_plot,
       width = 10, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# 10. Hidden Genes Analysis
# ==============================================================================

cat("\n\n========================================\n")
cat("'Hidden' Genes: Pathway-Supported but Not GWAS/MAGMA Significant\n")
cat("========================================\n")

hidden_genes <- gene_evidence %>%
  filter(
    hit_pathway == TRUE,
    hit_gwas == FALSE,
    hit_magma == FALSE
  ) %>%
  arrange(magma_p) %>%
  head(20)

if (nrow(hidden_genes) > 0) {
  cat("Top 20 'hidden' pathway genes (sorted by MAGMA p):\n")
  cat(strrep("-", 70), "\n")
  for (i in seq_len(nrow(hidden_genes))) {
    g <- hidden_genes[i, ]
    cat(
      sprintf("%2d. %s\n", i, g$GENE),
      sprintf("    MAGMA p: %s | GWAS min p: %s\n",
              format(g$magma_p, digits = 3, scientific = TRUE),
              ifelse(is.na(g$gwas_min_p), "NA",
                     format(g$gwas_min_p, digits = 3, scientific = TRUE))),
      sprintf("    Pathways: %s\n", g$pathways),
      "\n", sep = ""
    )
  }
} else {
  cat("No hidden genes found at current thresholds.\n")
}

# ==============================================================================
# 11. Core Genes Analysis
# ==============================================================================

cat("\n\n========================================\n")
cat("Core Genes: Present in Multiple Top Pathways\n")
cat("========================================\n")

core_genes <- gene_evidence %>%
  filter(n_top_pathways >= 2) %>%
  arrange(desc(n_top_pathways), magma_p)

if (nrow(core_genes) > 0) {
  cat("Genes in 2+ top pathways:", nrow(core_genes), "\n\n")
  for (i in seq_len(min(20, nrow(core_genes)))) {
    g <- core_genes[i, ]
    cat(
      sprintf("%2d. %s (in %d pathways)\n", i, g$GENE, g$n_top_pathways),
      sprintf("    MAGMA p: %s | GWAS min p: %s\n",
              format(g$magma_p, digits = 3, scientific = TRUE),
              ifelse(is.na(g$gwas_min_p), "NA",
                     format(g$gwas_min_p, digits = 3, scientific = TRUE))),
      sprintf("    Pathways: %s\n\n", g$pathways),
      sep = ""
    )
  }
} else {
  cat("No genes found in multiple top pathways.\n")
}

# ==============================================================================
# 12. Summary Visualization
# ==============================================================================

# Create summary bar chart with same categories as scatter plot
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
    sum(gene_evidence$hit_gwas & gene_evidence$hit_magma &
          gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & gene_evidence$hit_pathway &
          !gene_evidence$hit_magma, na.rm = TRUE),
    sum(gene_evidence$hit_magma & gene_evidence$hit_pathway &
          !gene_evidence$hit_gwas, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & gene_evidence$hit_magma &
          !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_gwas & !gene_evidence$hit_magma &
          !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_magma & !gene_evidence$hit_gwas &
          !gene_evidence$hit_pathway, na.rm = TRUE),
    sum(gene_evidence$hit_pathway & !gene_evidence$hit_gwas &
          !gene_evidence$hit_magma, na.rm = TRUE)
  )
)

summary_counts$Category <- factor(
  summary_counts$Category,
  levels = summary_counts$Category
)

summary_bar <- ggplot(summary_counts, aes(x = Category, y = Count, fill = Category)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = Count), vjust = -0.5, size = 4) +
  scale_fill_manual(values = c(
    "All 3" = "red",
    "GWAS + Pathway" = "hotpink",
    "MAGMA + Pathway" = "darkcyan",
    "GWAS + MAGMA" = "purple",
    "GWAS only" = "orange",
    "MAGMA only" = "navy",
    "Pathway only" = "forestgreen"
  )) +
  labs(x = "", y = "Number of Genes",
       title = "Candidate Gene Evidence Summary (Sorghum Stem Volume)") +
  theme_minimal() +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
    legend.position = "none"
  ) +
  ylim(0, max(summary_counts$Count) * 1.15)

ggsave(file.path(OUT_DIR, "candidate_genes_summary_sorghum_stem_vol.png"), summary_bar,
       width = 10, height = 8, dpi = 300, bg = "white")

# ==============================================================================
# 13. Final Summary
# ==============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("FINAL SUMMARY: Candidate Gene Analysis (Sorghum Stem Volume)\n")
cat(strrep("=", 70), "\n\n")

cat("Analysis Parameters:\n")
cat("  - SNP-to-gene window:", WINDOW_SIZE / 1000, "kb\n")
if (GWAS_MODE == "top_pct") {
  cat("  - GWAS: top", GWAS_TOP_PCT, "% of genes\n")
} else {
  cat("  - GWAS: p <", GWAS_P_THRESHOLD, "\n")
}
if (MAGMA_MODE == "top_pct") {
  cat("  - MAGMA: top", MAGMA_TOP_PCT, "% of genes\n")
} else {
  cat("  - MAGMA: FDR <", MAGMA_FDR_THRESHOLD, "\n")
}
if (PATHWAY_MODE == "fdr") {
  cat("  - Pathways: FDR <", PATHWAY_FDR_THRESHOLD,
      "(n =", n_pathways_used, ")\n")
} else {
  cat("  - Pathways: top", TOP_K_PATHWAYS, "\n")
}
cat("  - Gene universe (MAGMA-tested):", nrow(magma_gene), "\n\n")

cat("Evidence Layer Summary:\n")
cat(" ", gwas_label, ":",
    sum(gene_evidence$hit_gwas, na.rm = TRUE), "genes\n")
cat(" ", magma_label, ":",
    sum(gene_evidence$hit_magma, na.rm = TRUE), "genes\n")
cat(" ", pathway_label, ":",
    sum(gene_evidence$hit_pathway, na.rm = TRUE), "genes\n\n")

cat("Multi-Layer Support:\n")
cat("  Supported by 2+ layers:",
    sum(gene_evidence$support_layers >= 2, na.rm = TRUE), "genes\n")
cat("  Supported by all 3 layers:",
    sum(gene_evidence$support_layers >= 3, na.rm = TRUE), "genes\n\n")

hidden_count <- sum(gene_evidence$hit_pathway & !gene_evidence$hit_magma,
                    na.rm = TRUE)
cat("Key Finding:\n")
cat("  ", hidden_count, " genes are in significant pathways but NOT individually\n")
cat("  significant by", magma_label, ".\n\n")

cat("Pathway Gene Enrichment:\n")
cat("  Wilcoxon test p-value:", format(wilcox_result$p.value, digits = 4), "\n")

cat("\n========== DONE (Sorghum Stem Volume) ==========\n")
cat("Output files:\n")
cat("  -", file.path(OUT_DIR, "candidate_genes_top200_sorghum_stem_vol.csv"), "\n")
cat("  -", file.path(OUT_DIR, "candidate_genes_upset_sorghum_stem_vol.png"), "\n")
cat("  -", file.path(OUT_DIR, "pathway_gene_enrichment_density_sorghum_stem_vol.png"), "\n")
cat("  -", file.path(OUT_DIR, "gwas_vs_magma_scatter_sorghum_stem_vol.png"), "\n")
cat("  -", file.path(OUT_DIR, "candidate_genes_summary_sorghum_stem_vol.png"), "\n")
