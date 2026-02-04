# ==============================================================================
# Candidate Gene Analysis: Integrating GWAS, MAGMA, and Pathway Evidence
# Arabidopsis - Cold Temperature Phenotype
# ==============================================================================
#
# Three layers of evidence:
#   1. GWAS: SNP-level association (mapped to genes via 25kb window)
#      - Uses minimum SNP p-value per gene (most significant)
#      - Duplicate genes: keeps entry with lowest p-value
#
#   2. MAGMA: Gene-level association (LD-aware multi-SNP aggregation)
#      - Uses P_MULTI column (multi-SNP model p-value)
#      - FDR and Bonferroni corrections calculated here
#      - Duplicate genes: keeps entry with lowest p-value
#
#   3. Pathway: Set-level association (coordinated/diffuse signal detection)
#      - Genes from top K pathways by OMNIBUS p-value
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

GWAS_P_THRESHOLD   <- 1e-7     # GWAS SNP p-value threshold for gene hits
MAGMA_FDR_THRESHOLD <- 0.005    # MAGMA FDR (Benjamini-Hochberg) threshold
TOP_K_PATHWAYS      <- 10       # Number of top pathways to include
WINDOW_SIZE         <- 25000   # SNP-to-gene mapping window (bp)


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
    axis.title.x   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,
      color = "black"
    ),
    axis.line      = element_line(color = "black"),

    ## <-- CHANGED (was element_blank())
    panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 24)
  )

# ==============================================================================
# 1. Load Gene Location File (for SNP-to-gene mapping)
# ==============================================================================

gene_loc <- read.delim(
"/Users/nirwantandukar/Documents/Github/MAGCAT/inst/extdata/at.genes.loc",
stringsAsFactors = FALSE
)

head(gene_loc)
cat("Total genes in annotation:", nrow(gene_loc), "\n")

# Create extended gene regions
gene_loc <- gene_loc %>%
mutate(
  START_EXT = pmax(0, START - WINDOW_SIZE),
  STOP_EXT = STOP + WINDOW_SIZE
)

# ==============================================================================
# 2. Load GWAS Results and Map SNPs to Genes
# ==============================================================================

gwas_raw <- fread(
"/Users/nirwantandukar/Documents/Research/data/Arabidopsis/raw_gwas/Bio6_coldest_temp.csv",
header = TRUE
)

cat("Total SNPs in GWAS:", nrow(gwas_raw), "\n")

# Rename columns for clarity
names(gwas_raw)[names(gwas_raw) == "P-Value"] <- "P"
names(gwas_raw)[names(gwas_raw) == "Positions"] <- "POS"

# Check p-value distribution
cat("\nGWAS p-value summary:\n")
print(summary(gwas_raw$P))

# Count significant SNPs at different thresholds
cat("\nSNPs at different thresholds:\n")
cat("  p < 1e-7:", sum(gwas_raw$P < 1e-7, na.rm = TRUE), "\n")
cat("  p < 1e-6:", sum(gwas_raw$P < 1e-6, na.rm = TRUE), "\n")
cat("  p < 1e-5:", sum(gwas_raw$P < 1e-5, na.rm = TRUE), "\n")
cat("  p < 1e-4:", sum(gwas_raw$P < 1e-4, na.rm = TRUE), "\n")
cat("  p < 1e-3:", sum(gwas_raw$P < 1e-3, na.rm = TRUE), "\n")
cat("  p < 0.05:", sum(gwas_raw$P < 0.05, na.rm = TRUE), "\n")

# Map SNPs to genes using 25kb window
# For each SNP, find all genes where SNP falls within extended region

map_snps_to_genes <- function(gwas_df, gene_df, window = 25000) {
# This can be slow for large datasets, so we do it chromosome by chromosome

result_list <- list()

for (chr in unique(gwas_df$CHR)) {
  gwas_chr <- gwas_df[gwas_df$CHR == chr, ]
  gene_chr <- gene_df[gene_df$CHR == chr, ]

  if (nrow(gene_chr) == 0) next

  # For each SNP, find overlapping genes
  for (i in seq_len(nrow(gwas_chr))) {
    snp_pos <- gwas_chr$POS[i]
    snp_p <- gwas_chr$P[i]
    snp_id <- gwas_chr$SNP[i]

    # Find genes where SNP is within extended region
    overlaps <- gene_chr$GENE[
      snp_pos >= gene_chr$START_EXT & snp_pos <= gene_chr$STOP_EXT
    ]

    if (length(overlaps) > 0) {
      result_list[[length(result_list) + 1]] <- data.frame(
        SNP = snp_id,
        CHR = chr,
        POS = snp_pos,
        P = snp_p,
        GENE = overlaps,
        stringsAsFactors = FALSE
      )
    }
  }
}

do.call(rbind, result_list)
}

cat("\nMapping SNPs to genes (25kb window)... this may take a moment\n")

# For speed, let's use a vectorized approach with data.table
gwas_dt <- as.data.table(gwas_raw)
gene_dt <- as.data.table(gene_loc)

# Set keys for fast overlap joins
setkey(gene_dt, CHR, START_EXT, STOP_EXT)

# Use foverlaps for fast interval overlap
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

# Get minimum p-value per gene (best GWAS evidence for each gene)
# NOTE: If same gene appears multiple times (duplicate entries), we keep the
# entry with the LOWEST p-value (most significant association)
gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),  # Lowest p = most significant
    gwas_n_snps = n(),
    gwas_best_snp = SNP[which.min(P)],
    .groups = "drop"
  ) %>%
  # Remove any duplicate gene names, keeping the one with lowest p-value
  arrange(gwas_min_p) %>%
  filter(!duplicated(GENE))

cat("\nGenes with GWAS evidence:", nrow(gwas_gene), "\n")

# GWAS hits at different thresholds
cat("\nGenes with GWAS hit (min SNP p-value):\n")
cat("  p < 1e-7:", sum(gwas_gene$gwas_min_p < 1e-7), "\n")
cat("  p < 1e-6:", sum(gwas_gene$gwas_min_p < 1e-6), "\n")
cat("  p < 1e-5:", sum(gwas_gene$gwas_min_p < 1e-5), "\n")
cat("  p < 1e-4:", sum(gwas_gene$gwas_min_p < 1e-4), "\n")
cat("  p < 1e-3:", sum(gwas_gene$gwas_min_p < 1e-3), "\n")

# ==============================================================================
# 3. Load MAGMA Gene Results
# ==============================================================================

magma_files <- list.files(
path = "/Users/nirwantandukar/Documents/Research/results/Arabidopsis/MAGMA/AT_cold_by_chr",
pattern = "^AT_cold_chr.*\\.multi_snp_wise.genes\\.out$",
full.names = TRUE
)

cat("\nLoading MAGMA files:", length(magma_files), "\n")

magma_list <- lapply(magma_files, function(f) {
read.table(f, header = TRUE, stringsAsFactors = FALSE, comment.char = "#")
})

magma_all <- do.call(rbind, magma_list)

# Check columns available
cat("MAGMA columns:", paste(names(magma_all), collapse = ", "), "\n")

# Use P_MULTI as the primary MAGMA p-value
# P_MULTI = multi-SNP model p-value (preferred)
# P_SNPWISE_MEAN = mean of SNP-wise p-values
# P_SNPWISE_TOP1 = top SNP p-value
if (!"P_MULTI" %in% names(magma_all)) {
  stop("P_MULTI column not found in MAGMA output!")
}

# Rename P_MULTI to magma_p for clarity
magma_all$magma_p <- magma_all$P_MULTI

# Remove duplicates (keep best p-value per gene if same gene appears multiple times)
magma_gene <- magma_all %>%
  group_by(GENE) %>%
  slice_min(magma_p, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(GENE, CHR, START, STOP, NSNPS, ZSTAT, magma_p)

cat("Genes with MAGMA results:", nrow(magma_gene), "\n")

# Calculate multiple testing corrections (not in original MAGMA output)
# FDR = Benjamini-Hochberg False Discovery Rate
# Bonferroni = conservative family-wise error rate control
magma_gene <- magma_gene %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_bonf = p.adjust(magma_p, method = "bonferroni")
  ) %>%
  arrange(magma_p)

# MAGMA hits at different significance thresholds
cat("\nMAGMA significant genes (using P_MULTI):\n")
cat("  Bonferroni < 0.05:", sum(magma_gene$magma_bonf < 0.05, na.rm = TRUE), "\n")
cat("  FDR < 0.05:", sum(magma_gene$magma_fdr < 0.05, na.rm = TRUE), "\n")
cat("  FDR < 0.10:", sum(magma_gene$magma_fdr < 0.10, na.rm = TRUE), "\n")
cat("  Nominal p < 0.05:", sum(magma_gene$magma_p < 0.05, na.rm = TRUE), "\n")
cat("  Top 100 genes (for ranking)\n")

# Show top 10 MAGMA genes
cat("\nTop 10 MAGMA genes by P_MULTI:\n")
head(magma_gene, 10) %>%
  select(GENE, NSNPS, ZSTAT, magma_p, magma_fdr, magma_bonf) %>%
  print()

# ==============================================================================
# 4. Load Pathway Results
# ==============================================================================

omni_results <- readRDS(
"/Users/nirwantandukar/Documents/Github/MAGCAT/omni_results_arabidopsis_ACAT_B10000.RDS"
)

cat("\nTotal pathways:", nrow(omni_results), "\n")

# Get omnibus p-value column
omni_col <- "omni_p_final"
if (!omni_col %in% names(omni_results)) {
omni_col <- "omni_p_mvn"
if (!omni_col %in% names(omni_results)) {
  omni_col <- "omni_p_analytic"
}
}
cat("Using pathway p-value column:", omni_col, "\n")

# Check pathway hits
cat("\nPathway significant counts:\n")
bonf_thresh <- 0.05 / nrow(omni_results)
cat("  Bonferroni (", format(bonf_thresh, digits = 3), "):",
  sum(omni_results[[omni_col]] < bonf_thresh, na.rm = TRUE), "\n")
cat("  FDR < 0.10:",
  sum(p.adjust(omni_results[[omni_col]], "BH") < 0.10, na.rm = TRUE), "\n")
cat("  Nominal < 0.05:",
  sum(omni_results[[omni_col]] < 0.05, na.rm = TRUE), "\n")

# Get top K pathways
top_pathways <- omni_results %>%
arrange(!!sym(omni_col)) %>%
slice_head(n = TOP_K_PATHWAYS)

cat("\nTop", TOP_K_PATHWAYS, "pathways:\n")
for (i in seq_len(nrow(top_pathways))) {
cat("  ", i, ". ", top_pathways$pathway_id[i], " (",
    top_pathways$pathway_name[i], ")\n",
    "     p = ", format(top_pathways[[omni_col]][i], digits = 3), "\n",
    sep = "")
}

# Extract genes from top pathways
# The genes_used column contains semicolon-separated gene lists
pathway_genes_list <- strsplit(top_pathways$genes_used, ";")
names(pathway_genes_list) <- top_pathways$pathway_id

# Create pathway_genes data frame
pathway_genes <- data.frame(
pathway_id = rep(top_pathways$pathway_id, lengths(pathway_genes_list)),
pathway_name = rep(top_pathways$pathway_name, lengths(pathway_genes_list)),
pathway_p = rep(top_pathways[[omni_col]], lengths(pathway_genes_list)),
GENE = unlist(pathway_genes_list),
stringsAsFactors = FALSE
)

# Trim whitespace from gene names
pathway_genes$GENE <- trimws(pathway_genes$GENE)

cat("\nGenes in top", TOP_K_PATHWAYS, "pathways:", nrow(pathway_genes), "\n")
cat("Unique genes in top pathways:",
  length(unique(pathway_genes$GENE)), "\n")

# Summarize pathway support per gene
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
cat("Integrating Evidence Across Layers\n")
cat("========================================\n")

# Define the gene universe as MAGMA-tested genes
gene_universe <- unique(magma_gene$GENE)
cat("\nGene universe (MAGMA-tested genes):", length(gene_universe), "\n")

# Merge all evidence
gene_evidence <- magma_gene %>%
left_join(gwas_gene, by = "GENE") %>%
left_join(pathway_gene_support, by = "GENE") %>%
mutate(
  # Define hit criteria using top-level parameters
  hit_gwas    = !is.na(gwas_min_p) & gwas_min_p < GWAS_P_THRESHOLD,
  hit_magma   = !is.na(magma_fdr) & magma_fdr < MAGMA_FDR_THRESHOLD,
  hit_pathway = !is.na(n_top_pathways) & n_top_pathways >= 1,

  # Support layers count
  support_layers = hit_gwas + hit_magma + hit_pathway,

  # Composite score
  score = (hit_gwas * 1) + (hit_magma * 1) + (hit_pathway * 1) +
    0.2 * ifelse(!is.na(magma_p), -log10(magma_p), 0) +
    0.1 * ifelse(!is.na(gwas_min_p), -log10(gwas_min_p), 0) +
    0.1 * ifelse(!is.na(best_pathway_p), -log10(best_pathway_p), 0)
) %>%
arrange(desc(score))

# ==============================================================================
# 6. Summary Statistics
# ==============================================================================

cat("\n\nHit counts per layer:\n")
cat("  GWAS (p <", GWAS_P_THRESHOLD, "):",
  sum(gene_evidence$hit_gwas, na.rm = TRUE), "\n")
cat("  MAGMA (FDR <", MAGMA_FDR_THRESHOLD, "):",
  sum(gene_evidence$hit_magma, na.rm = TRUE), "\n")
cat("  Pathway (top", TOP_K_PATHWAYS, "):",
  sum(gene_evidence$hit_pathway, na.rm = TRUE), "\n")

cat("\nMulti-layer support:\n")
for (n in 0:3) {
count <- sum(gene_evidence$support_layers == n, na.rm = TRUE)
cat("  ", n, " layers:", count, "\n")
}

# ==============================================================================
# 7. UpSet Plot: Overlap of Evidence Layers
# ==============================================================================

# Create binary matrix for UpSet
gwas_label   <- paste0("GWAS (p<", GWAS_P_THRESHOLD, ")")
magma_label  <- paste0("MAGMA (FDR<", MAGMA_FDR_THRESHOLD, ")")
pathway_label <- paste0("Pathway (top ", TOP_K_PATHWAYS, ")")

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

# Filter to genes with at least one hit
upset_df <- upset_df[rowSums(upset_df) > 0, ]

cat("\nGenes with at least one layer of evidence:", nrow(upset_df), "\n")

png("candidate_genes_upset_arabidopsis.png",
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

cat("UpSet plot saved to candidate_genes_upset.pdf\n")

# ==============================================================================
# 8. Top Candidate Genes Table
# ==============================================================================

# Get top candidates (genes with score > 0 or multi-layer support)
top_candidates <- gene_evidence %>%
filter(score > 1 | support_layers >= 2) %>%
arrange(desc(score)) %>%
head(50)

cat("\n\n========================================\n")
cat("TOP CANDIDATE GENES\n")
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

# Save full table
write.csv(
gene_evidence %>%
  select(GENE, magma_p, magma_fdr, gwas_min_p, gwas_n_snps,
         n_top_pathways, best_pathway_p, pathways,
         hit_gwas, hit_magma, hit_pathway,
         support_layers, score) %>%
  head(200),
"candidate_genes_top200_arabidopsis.csv",
row.names = FALSE
)

cat("\nFull table saved to candidate_genes_top200_arabidopsis.csv\n")

# ==============================================================================
# 9. Pathway Gene Enrichment Analysis
# ==============================================================================

cat("\n\n========================================\n")
cat("Pathway Gene Enrichment Analysis\n")
cat("========================================\n")
cat("Are genes in significant pathways enriched for lower MAGMA p-values?\n\n")

# Compare MAGMA p-values: pathway genes vs non-pathway genes
gene_evidence$in_pathway <- !is.na(gene_evidence$n_top_pathways)

pathway_genes_magma <- gene_evidence$magma_p[gene_evidence$in_pathway]
nonpathway_genes_magma <- gene_evidence$magma_p[!gene_evidence$in_pathway]

cat("Genes in top pathways:", length(pathway_genes_magma), "\n")
cat("Genes not in top pathways:", length(nonpathway_genes_magma), "\n")

# Wilcoxon test
wilcox_result <- wilcox.test(
pathway_genes_magma,
nonpathway_genes_magma,
alternative = "less"  # pathway genes have LOWER p-values
)

cat("\nWilcoxon rank-sum test (pathway genes have lower MAGMA p?):\n")
cat("  W =", wilcox_result$statistic, "\n")
cat("  p-value =", format(wilcox_result$p.value, digits = 4), "\n")

# Summary statistics
cat("\nMedian -log10(MAGMA p):\n")
cat("  Pathway genes:", round(median(-log10(pathway_genes_magma), na.rm = TRUE), 3), "\n")
cat("  Non-pathway genes:", round(median(-log10(nonpathway_genes_magma), na.rm = TRUE), 3), "\n")

# Density plot
density_df <- gene_evidence %>%
mutate(Group = ifelse(in_pathway, "In Top Pathways", "Not in Top Pathways"))

density_plot <- ggplot(density_df, aes(x = -log10(magma_p), fill = Group)) +
geom_density(alpha = 0.5) +
labs(
  title = "MAGMA P-value Distribution: Pathway vs Non-Pathway Genes",
  subtitle = paste0("Wilcoxon p = ", format(wilcox_result$p.value, digits = 3)),
  x = expression(-log[10](MAGMA~p)),
  y = "Density"
) +
scale_fill_manual(values = c("forestgreen", "gray50")) +
theme_minimal() +
theme(legend.position = "bottom") + plot_theme

quartz()
print(density_plot)

#ggsave("pathway_gene_enrichment_density.pdf", density_plot, width = 10, height = 6)
ggsave("pathway_gene_enrichment_density.png", density_plot, width = 8, height = 8,
     dpi = 300, bg = "white")

# ==============================================================================
# 10. Scatter Plot: GWAS vs MAGMA with Pathway Annotation
# ==============================================================================

# Filter to genes with both GWAS and MAGMA evidence
scatter_df <- gene_evidence %>%
filter(!is.na(gwas_min_p) & !is.na(magma_p)) %>%
mutate(
  highlight = case_when(
    hit_gwas & hit_magma & hit_pathway ~ "All 3 layers",
    (hit_gwas | hit_magma) & hit_pathway ~ "2 layers + Pathway",
    hit_gwas & hit_magma ~ "GWAS + MAGMA only",
    hit_pathway ~ "Pathway only",
    TRUE ~ "None significant"
  ),
  highlight = factor(highlight, levels = c(
    "All 3 layers", "2 layers + Pathway", "GWAS + MAGMA only",
    "Pathway only", "None significant"
  ))
)

scatter_plot <- ggplot(
scatter_df,
aes(x = -log10(gwas_min_p), y = -log10(magma_p), color = highlight)
) +
geom_point(alpha = 0.6, size = 2) +
geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "gray50") +
geom_vline(xintercept = 4, linetype = "dashed", color = "gray50") +
scale_color_manual(values = c(
  "All 3 layers" = "red",
  "2 layers + Pathway" = "orange",
  "GWAS + MAGMA only" = "purple",
  "Pathway only" = "forestgreen",
  "None significant" = "gray70"
)) +
labs(
  title = "GWAS vs MAGMA Evidence with Pathway Annotation",
  x = expression(-log[10](GWAS~min~p)),
  y = expression(-log[10](MAGMA~p)),
  color = "Evidence"
) +
theme_minimal() +
theme(legend.position = "right") + plot_theme +
  guides(
    color = guide_legend(
      nrow = 3,       # 2 rows
      ncol = 2,       # 3 columns
      byrow = TRUE,   # fill across rows first (usually what you want)
      override.aes = list(size = 4, alpha = 1)  # make legend dots visible
    )
  )

quartz()
print(scatter_plot)
#ggsave("gwas_vs_magma_scatter.pdf", scatter_plot, width = 10, height = 8)
ggsave("gwas_vs_magma_scatter_arabidopsis.png", scatter_plot, width = 8, height = 8,
     dpi = 300, bg = "white")

# ==============================================================================
# 11. "Hidden" Genes: Pathway-Supported but Missed by GWAS/MAGMA
# ==============================================================================

cat("\n\n========================================\n")
cat("'Hidden' Genes: Pathway-Supported but Not GWAS/MAGMA Significant\n")
cat("========================================\n")
cat("These genes are in significant pathways but don't reach individual\n")
cat("significance - demonstrating pathway analysis captures diffuse signal.\n\n")

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
# 12. Core Genes Analysis (Genes in Multiple Top Pathways)
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
# 13. Summary Visualization
# ==============================================================================

# Create summary bar chart
summary_counts <- data.frame(
Category = c(
  gwas_label,
  magma_label,
  pathway_label,
  "2+ layers",
  "All 3 layers"
),
Count = c(
  sum(gene_evidence$hit_gwas, na.rm = TRUE),
  sum(gene_evidence$hit_magma, na.rm = TRUE),
  sum(gene_evidence$hit_pathway, na.rm = TRUE),
  sum(gene_evidence$support_layers >= 2, na.rm = TRUE),
  sum(gene_evidence$support_layers >= 3, na.rm = TRUE)
),
Type = c(
  "Single Layer", "Single Layer", "Single Layer",
  "Multi-Layer", "Multi-Layer"
)
)

summary_counts$Category <- factor(
summary_counts$Category,
levels = summary_counts$Category
)

summary_bar <- ggplot(summary_counts, aes(x = Category, y = Count, fill = Type)) +
geom_col(alpha = 0.8) +
geom_text(aes(label = Count), vjust = -0.5, size = 4) +
scale_fill_manual(values = c("Single Layer" = "steelblue",
                             "Multi-Layer" = "darkred")) +
labs(
  #title = "Candidate Gene Counts by Evidence Layer",
  #subtitle = "Multi-layer support identifies high-confidence candidates",
  x = "",
  y = "Number of Genes"
) +
theme_minimal() +
theme(
  axis.text.x = element_text(angle = 45, hjust = 1),
  legend.position = "bottom"
) +
ylim(0, max(summary_counts$Count) * 1.15) + plot_theme

summary_bar <- summary_bar +
  theme(
    axis.text.x = element_text(
      angle = 35,      # slant (try 25–45)
      hjust = 1,
      vjust = 1
    )
  ) + theme(legend.title = element_blank())

quartz()
print(summary_bar)

#ggsave("candidate_genes_summary.pdf", summary_bar, width = 10, height = 6)
ggsave("candidate_genes_summary_arabidopsis.png", summary_bar, width = 8, height = 8,
     dpi = 300, bg = "white")

# ==============================================================================
# 14. Final Summary
# ==============================================================================

cat("\n\n")
cat(strrep("=", 70), "\n")
cat("FINAL SUMMARY: Candidate Gene Analysis\n")
cat(strrep("=", 70), "\n\n")

cat("Analysis Parameters:\n")
cat("  - SNP-to-gene window:", WINDOW_SIZE / 1000, "kb\n")
cat("  - GWAS p threshold:", GWAS_P_THRESHOLD, "\n")
cat("  - MAGMA FDR threshold:", MAGMA_FDR_THRESHOLD, "\n")
cat("  - Top pathways used:", TOP_K_PATHWAYS, "\n")
cat("  - Gene universe (MAGMA-tested):", nrow(magma_gene), "\n\n")

cat("Evidence Layer Summary:\n")
cat("  GWAS (p <", GWAS_P_THRESHOLD, "):",
  sum(gene_evidence$hit_gwas, na.rm = TRUE), "genes\n")
cat("  MAGMA (FDR <", MAGMA_FDR_THRESHOLD, "):",
  sum(gene_evidence$hit_magma, na.rm = TRUE), "genes\n")
cat("  Pathway (top", TOP_K_PATHWAYS, "):",
  sum(gene_evidence$hit_pathway, na.rm = TRUE), "genes\n\n")

cat("Multi-Layer Support:\n")
cat("  Supported by 2+ layers:",
  sum(gene_evidence$support_layers >= 2, na.rm = TRUE), "genes\n")
cat("  Supported by all 3 layers:",
  sum(gene_evidence$support_layers >= 3, na.rm = TRUE), "genes\n\n")

cat("Key Finding:\n")
hidden_count <- sum(gene_evidence$hit_pathway & !gene_evidence$hit_magma,
                  na.rm = TRUE)
cat("  ", hidden_count, " genes are in significant pathways but NOT individually\n")
cat("  significant by MAGMA (FDR <", MAGMA_FDR_THRESHOLD, ").\n")
cat("  This demonstrates that pathway analysis captures coordinated/diffuse\n")
cat("  signals missed by gene-level tests.\n\n")

cat("Pathway Gene Enrichment:\n")
cat("  Wilcoxon test p-value:", format(wilcox_result$p.value, digits = 4), "\n")
if (wilcox_result$p.value < 0.05) {
cat("  ** Pathway genes have significantly lower MAGMA p-values,\n")
cat("     validating that pathway hits are enriched for genetic signal. **\n")
}

cat("\n========== DONE ==========\n")
cat("Output files:\n")
cat("  - candidate_genes_top200.csv\n")
cat("  - candidate_genes_upset.pdf\n")
cat("  - pathway_gene_enrichment_density.pdf/png\n")
cat("  - gwas_vs_magma_scatter.pdf/png\n")
cat("  - candidate_genes_summary.pdf/png\n")

