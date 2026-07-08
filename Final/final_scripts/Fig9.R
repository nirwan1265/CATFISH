#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(UpSetR)
})

default_repo <- "/Users/nirwantandukar/Documents/Github/MAGCAT"
default_results_dir <- file.path(default_repo, "Final", "final_main_fig", "SAP_stem_volume_realdata")

GWAS_MODE <- "top_pct"
GWAS_TOP_PCT <- 1
GWAS_P_THRESHOLD <- 1e-5

MAGMA_MODE <- "top_pct"
MAGMA_TOP_PCT <- 5
MAGMA_FDR_THRESHOLD <- 0.3

PATHWAY_MODE <- "fdr"
PATHWAY_FDR_THRESHOLD <- 0.05
TOP_K_PATHWAYS <- 20
WINDOW_SIZE <- 25000

RESULTS_DIR <- Sys.getenv("RESULTS_DIR", default_results_dir)
GWAS_FILE <- Sys.getenv(
  "GWAS_FILE",
  "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt"
)
GENE_LOC_FILE <- Sys.getenv("GENE_LOC_FILE", file.path(default_repo, "inst/extdata/sorghum.genes.loc"))

pick_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else NA_character_
}

catfish_csv <- pick_existing_file(c(
  file.path(RESULTS_DIR, "sorghum_stem_vol_catfish_results_MVN.csv"),
  file.path(RESULTS_DIR, "sorghum_stem_vol_CATFISH_B10000_GPD.csv"),
  file.path(RESULTS_DIR, "sorghum_stem_vol_CATFISH_B1000000_GPD.csv")
))
genes_file <- pick_existing_file(c(
  file.path(RESULTS_DIR, "sorghum_stem_vol_genes_combined.txt"),
  file.path(RESULTS_DIR, "Stem_volume_combined_genes.tsv")
))

if (!file.exists(GWAS_FILE)) stop("Missing GWAS file: ", GWAS_FILE, call. = FALSE)
if (!file.exists(GENE_LOC_FILE)) stop("Missing gene location file: ", GENE_LOC_FILE, call. = FALSE)
if (!nzchar(catfish_csv) || !file.exists(catfish_csv)) stop("Missing CATFISH CSV in ", RESULTS_DIR, call. = FALSE)
if (!nzchar(genes_file) || !file.exists(genes_file)) stop("Missing combined genes file in ", RESULTS_DIR, call. = FALSE)

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

plot_theme <- theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 14, color = "black"),
    axis.line = element_line(color = "black"),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.35),
    panel.grid.minor = element_line(color = "grey95", linewidth = 0.25),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

method_cols <- c(
  ACAT = "acat_p",
  Fisher = "fisher_p",
  TFisher = "tfisher_p_analytic",
  minP = "minp_p_analytic",
  Stouffer = "stouffer_p_analytic",
  Omnibus = "omni_p_final"
)

method_colors <- c(
  ACAT = "#D95F02",
  Fisher = "#CC79A7",
  TFisher = "#6C63FF",
  minP = "#E6AB02",
  Stouffer = "#1B9E77",
  Omnibus = "#C9332B"
)

omni_results <- read.csv(catfish_csv, stringsAsFactors = FALSE)
genes <- fread(genes_file, header = TRUE)
gene_loc <- read.delim(GENE_LOC_FILE, stringsAsFactors = FALSE)
gwas_raw <- fread(GWAS_FILE, header = TRUE)

names(gwas_raw)[names(gwas_raw) == "p_wald"] <- "P"
names(gwas_raw)[names(gwas_raw) == "ps"] <- "POS"
names(gwas_raw)[names(gwas_raw) == "chr"] <- "CHR"
names(gwas_raw)[names(gwas_raw) == "rs"] <- "SNP_ID"

gene_loc <- gene_loc %>%
  mutate(
    START_EXT = pmax(0, START - WINDOW_SIZE),
    STOP_EXT = STOP + WINDOW_SIZE
  )

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

gwas_gene <- snp_gene_map %>%
  as.data.frame() %>%
  group_by(GENE) %>%
  summarise(
    gwas_min_p = min(P, na.rm = TRUE),
    gwas_n_snps = n(),
    .groups = "drop"
  ) %>%
  arrange(gwas_min_p) %>%
  mutate(gwas_rank = rank(gwas_min_p, ties.method = "min"))

genes <- as.data.frame(genes)
if (!"P" %in% names(genes) && "P_MULTI" %in% names(genes)) {
  genes$P <- genes$P_MULTI
}
genes <- genes %>%
  rename(magma_p = P) %>%
  mutate(
    magma_fdr = p.adjust(magma_p, method = "BH"),
    magma_bonf = p.adjust(magma_p, method = "bonferroni")
  ) %>%
  arrange(magma_p)

if (GWAS_MODE == "top_pct") {
  gwas_n_top <- ceiling(nrow(gwas_gene) * GWAS_TOP_PCT / 100)
  gwas_eff_threshold <- sort(gwas_gene$gwas_min_p)[gwas_n_top]
} else {
  gwas_eff_threshold <- GWAS_P_THRESHOLD
}

if (MAGMA_MODE == "top_pct") {
  magma_n_top <- ceiling(nrow(genes) * MAGMA_TOP_PCT / 100)
  magma_eff_threshold <- sort(genes$magma_p)[magma_n_top]
} else {
  magma_eff_threshold <- NULL
}

omni_col <- if ("omni_p_final" %in% names(omni_results)) {
  "omni_p_final"
} else if ("omni_p_mvn" %in% names(omni_results)) {
  "omni_p_mvn"
} else {
  "omni_p_analytic"
}

if (PATHWAY_MODE == "fdr") {
  omni_results$pathway_fdr <- p.adjust(omni_results[[omni_col]], "BH")
  top_pathways <- omni_results %>%
    filter(pathway_fdr < PATHWAY_FDR_THRESHOLD) %>%
    arrange(.data[[omni_col]])
} else {
  top_pathways <- omni_results %>%
    arrange(.data[[omni_col]]) %>%
    slice_head(n = TOP_K_PATHWAYS)
}

pathway_genes_list <- strsplit(top_pathways$genes_used, ";")
pathway_genes <- data.frame(
  pathway_id = rep(top_pathways$pathway_id, lengths(pathway_genes_list)),
  pathway_name = rep(top_pathways$pathway_name, lengths(pathway_genes_list)),
  pathway_p = rep(top_pathways[[omni_col]], lengths(pathway_genes_list)),
  GENE = trimws(unlist(pathway_genes_list)),
  stringsAsFactors = FALSE
)

pathway_gene_support <- pathway_genes %>%
  group_by(GENE) %>%
  summarise(
    n_top_pathways = n_distinct(pathway_id),
    best_pathway_p = min(pathway_p, na.rm = TRUE),
    pathways = paste(pathway_id, collapse = "; "),
    .groups = "drop"
  )

gene_evidence <- genes %>%
  left_join(gwas_gene, by = "GENE") %>%
  left_join(pathway_gene_support, by = "GENE") %>%
  mutate(
    hit_gwas = !is.na(gwas_min_p) & gwas_min_p <= gwas_eff_threshold,
    hit_magma = if (MAGMA_MODE == "top_pct") {
      !is.na(magma_p) & magma_p <= magma_eff_threshold
    } else {
      !is.na(magma_fdr) & magma_fdr < MAGMA_FDR_THRESHOLD
    },
    hit_pathway = !is.na(n_top_pathways) & n_top_pathways >= 1
  )

bonf_threshold <- 0.05 / nrow(omni_results)

# Panel A: UpSet
sig_sets <- lapply(method_cols, function(col) {
  omni_results$pathway_id[omni_results[[col]] < bonf_threshold]
})
names(sig_sets) <- names(method_cols)

all_pathways <- unique(unlist(sig_sets))
upset_df <- as.data.frame(sapply(names(sig_sets), function(s) {
  as.integer(all_pathways %in% sig_sets[[s]])
}))
rownames(upset_df) <- all_pathways

upset_png <- file.path(RESULTS_DIR, "tmp_sap_upset.png")
png(upset_png, width = 2400, height = 1800, res = 250)
upset(
  upset_df,
  sets = c("Omnibus", "ACAT", "Fisher", "TFisher", "minP", "Stouffer"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  mainbar.y.label = "Intersection Size",
  sets.x.label = "Set Size",
  text.scale = c(2.1, 2.0, 2.0, 1.8, 1.8, 1.8),
  point.size = 5,
  line.size = 1.6,
  mb.ratio = c(0.6, 0.4),
  sets.bar.color = c("#B22222", rep("#4C78A8", 5)),
  main.bar.color = "gray25",
  set_size.angles = 0
)
dev.off()

panel_a <- ggdraw() +
  draw_image(upset_png) +
  draw_label("A", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

# Panel B: rank curves
rank_long <- bind_rows(lapply(names(method_cols), function(method) {
  col <- method_cols[[method]]
  pvals <- omni_results[[col]]
  pvals <- pvals[is.finite(pvals) & !is.na(pvals) & pvals > 0 & pvals <= 1]
  data.frame(
    Method = method,
    Rank = seq_along(pvals),
    neglog10p = -log10(sort(pmax(pvals, 1e-15))),
    stringsAsFactors = FALSE
  )
}))

panel_b <- ggplot(rank_long, aes(Rank, neglog10p, color = Method)) +
  geom_line(linewidth = 1.1) +
  scale_color_manual(values = method_colors) +
  labs(
    x = "Pathway Rank",
    y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 10),
    plot.margin = margin(10, 10, 10, 10)
  )
panel_b <- ggdraw(panel_b) +
  draw_label("B", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

# Panel C: combined matrix
method_names <- names(method_cols)
n_methods <- length(method_names)
jaccard_mat <- matrix(NA_real_, n_methods, n_methods, dimnames = list(method_names, method_names))
cor_mat <- matrix(NA_real_, n_methods, n_methods, dimnames = list(method_names, method_names))

for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    set_i <- sig_sets[[method_names[i]]]
    set_j <- sig_sets[[method_names[j]]]
    inter <- length(intersect(set_i, set_j))
    uni <- length(union(set_i, set_j))
    jaccard_mat[i, j] <- if (uni > 0) inter / uni else 0

    p_i <- -log10(pmax(omni_results[[method_cols[[method_names[i]]]]], 1e-15))
    p_j <- -log10(pmax(omni_results[[method_cols[[method_names[j]]]]], 1e-15))
    valid <- is.finite(p_i) & is.finite(p_j) & !is.na(p_i) & !is.na(p_j)
    cor_mat[i, j] <- cor(p_i[valid], p_j[valid], method = "spearman")
  }
}

combined_mat <- matrix(NA_real_, n_methods, n_methods, dimnames = list(method_names, method_names))
for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    if (i < j) {
      combined_mat[i, j] <- cor_mat[i, j]
    } else if (i > j) {
      combined_mat[i, j] <- jaccard_mat[i, j]
    } else {
      combined_mat[i, j] <- 1
    }
  }
}

combined_long <- as.data.frame(as.table(combined_mat))
names(combined_long) <- c("Method1", "Method2", "Value")
combined_long$Method1 <- factor(combined_long$Method1, levels = method_names)
combined_long$Method2 <- factor(combined_long$Method2, levels = method_names)
combined_long$Label <- sprintf("%.2f", combined_long$Value)

panel_c_base <- ggplot(combined_long, aes(Method2, Method1, fill = Value)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = Label), size = 3.3, fontface = "bold") +
  scale_fill_gradient2(
    low = "#4575B4",
    mid = "white",
    high = "#B30000",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Value"
  ) +
  labs(
    x = "",
    y = "",
    caption = "Upper triangle: Spearman correlation of -log10(p) | Lower triangle: Jaccard similarity of Bonferroni-significant sets"
  ) +
  plot_theme +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    plot.caption = element_text(size = 9, hjust = 0.5),
    panel.grid = element_blank()
  ) +
  coord_fixed()
panel_c <- ggdraw(panel_c_base) +
  draw_label("C", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

# Panel D: GWAS vs MAGMA scatter
if (MAGMA_MODE == "top_pct") {
  magma_p_threshold <- magma_eff_threshold
} else {
  magma_p_threshold <- max(gene_evidence$magma_p[gene_evidence$magma_fdr < MAGMA_FDR_THRESHOLD], na.rm = TRUE)
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
    point_size = ifelse(highlight == "None", 1.6, 2.8)
  )

scatter_colors <- c(
  "All 3" = "#FF3B3B",
  "GWAS + Pathway" = "#16A3A6",
  "MAGMA + Pathway" = "#2E8B8B",
  "GWAS + MAGMA" = "#9C5DFF",
  "GWAS only" = "#3E5DB3",
  "MAGMA only" = "#3554B1",
  "Pathway only" = "#4CAF50",
  "None" = "#CFCFCF"
)

panel_d_base <- ggplot(
  scatter_df,
  aes(x = -log10(gwas_min_p), y = -log10(magma_p), color = highlight, size = point_size)
) +
  geom_point(alpha = 0.65) +
  scale_size_identity() +
  geom_hline(yintercept = -log10(magma_p_threshold), linetype = "dashed", color = "gray60") +
  geom_vline(xintercept = -log10(gwas_eff_threshold), linetype = "dashed", color = "gray60") +
  scale_color_manual(values = scatter_colors) +
  labs(
    x = expression(-log[10](GWAS~p)),
    y = expression(-log[10](MAGMA~p))
  ) +
  plot_theme +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  ) +
  guides(color = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.5, alpha = 1)))
panel_d <- ggdraw(panel_d_base) +
  draw_label("D", x = 0.01, y = 0.99, hjust = 0, vjust = 1, fontface = "bold", size = 20)

top_row <- plot_grid(panel_a, panel_b, ncol = 2, rel_widths = c(1.08, 1.0))
bottom_row <- plot_grid(panel_c, panel_d, ncol = 2, rel_widths = c(1.0, 1.0))
full_fig <- plot_grid(top_row, bottom_row, nrow = 2, rel_heights = c(1.0, 1.0))

png_file <- file.path(RESULTS_DIR, "Figure_SAP_Stem_volume_Method_Comparison.png")
pdf_file <- file.path(RESULTS_DIR, "Figure_SAP_Stem_volume_Method_Comparison.pdf")

ggsave(png_file, full_fig, width = 15, height = 12, dpi = 300, bg = "white")
ggsave(pdf_file, full_fig, width = 15, height = 12, bg = "white")

cat("Saved figure:\n")
cat("  ", png_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
