#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

default_repo <- Sys.getenv("CATFISH_REPO", "/Users/nirwantandukar/Documents/Github/MAGCAT")

get_env_int <- function(name, default) {
  x <- suppressWarnings(as.integer(Sys.getenv(name, "")))
  if (!is.finite(x) || is.na(x)) default else x
}

pick_existing_file <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit)) hit[[1]] else NA_character_
}

pick_final_p_col <- function(x) {
  picks <- c("omni_p_final", "omni_p_mvn", "omni_p_global", "omni_p_analytic")
  hit <- picks[picks %in% names(x)]
  if (length(hit)) hit[[1]] else NA_character_
}

TAU_OPTION <- Sys.getenv("TAU_OPTION", "paper")
TAU_LABEL <- Sys.getenv("TAU_LABEL", if (identical(TAU_OPTION, "paper")) "paper_tau_false" else TAU_OPTION)
B_PERM <- get_env_int("B_PERM", 10000L)
PERM_MODE <- Sys.getenv("PERM_MODE", "mvn")
TAIL_MODE <- Sys.getenv("TAIL_MODE", "hybrid_gpd")
tail_label <- if (identical(TAIL_MODE, "hybrid_gpd")) "GPD" else "empirical"

MAGMA_DIR <- Sys.getenv(
  "MAGMA_DIR",
  "/Users/nirwantandukar/Documents/Research/results/CATFISH/MAGMA/SAP_Stem_volume"
)
RESULTS_DIR <- Sys.getenv(
  "RESULTS_DIR",
  file.path(MAGMA_DIR, paste0("CATFISH_permutation_B", B_PERM, "_", PERM_MODE, "_", tail_label, "_", TAU_LABEL))
)
GWAS_FILE <- Sys.getenv(
  "GWAS_FILE",
  "/Users/nirwantandukar/Documents/Research/results/GWAS/SAP/Stem_diameter/Stem_volume_mod_sub_stem_volume_SAP_bialleles_MAF_0.05_11.assoc.txt"
)
OUT_PREFIX <- Sys.getenv("OUT_PREFIX", "Stem_volume")
TOP_N_PATHWAYS <- get_env_int("TOP_N_PATHWAYS", 18L)

combined_gene_file <- pick_existing_file(c(
  file.path(RESULTS_DIR, "sorghum_stem_vol_genes_combined.txt"),
  file.path(RESULTS_DIR, paste0(OUT_PREFIX, "_combined_genes.tsv"))
))

catfish_csv <- Sys.getenv("CATFISH_CSV", "")
if (!nzchar(catfish_csv)) {
  legacy_candidates <- c(
    file.path(RESULTS_DIR, "sorghum_stem_vol_catfish_results_MVN.csv"),
    file.path(RESULTS_DIR, paste0("sorghum_stem_vol_CATFISH_B", B_PERM, "_", tail_label, ".csv"))
  )
  catfish_csv <- pick_existing_file(legacy_candidates)
}
if (!nzchar(catfish_csv)) {
  cand <- list.files(
    RESULTS_DIR,
    pattern = paste0("^", OUT_PREFIX, "_CATFISH_.*\\.csv$"),
    full.names = TRUE
  )
  if (!length(cand)) {
    cand <- list.files(
      RESULTS_DIR,
      pattern = "^sorghum_stem_vol_.*\\.csv$",
      full.names = TRUE
    )
  }
  if (length(cand)) {
    catfish_csv <- cand[[1]]
  }
}

if (!file.exists(GWAS_FILE)) {
  stop("Missing GWAS_FILE: ", GWAS_FILE, call. = FALSE)
}
if (!nzchar(combined_gene_file) || !file.exists(combined_gene_file)) {
  stop("Missing combined MAGMA gene file in RESULTS_DIR: ", RESULTS_DIR, call. = FALSE)
}
if (!nzchar(catfish_csv) || !file.exists(catfish_csv)) {
  stop("Missing CATFISH CSV in RESULTS_DIR: ", RESULTS_DIR, call. = FALSE)
}

dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

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

cat("Using RESULTS_DIR: ", RESULTS_DIR, "\n", sep = "")
cat("Using GWAS_FILE: ", GWAS_FILE, "\n", sep = "")
cat("Using combined_gene_file: ", combined_gene_file, "\n", sep = "")
cat("Using catfish_csv: ", catfish_csv, "\n\n", sep = "")

gwas_raw <- fread(GWAS_FILE, header = TRUE)
names(gwas_raw)[names(gwas_raw) == "p_wald"] <- "P"
names(gwas_raw)[names(gwas_raw) == "ps"] <- "POS"
names(gwas_raw)[names(gwas_raw) == "chr"] <- "CHR"

gwas_raw <- gwas_raw %>%
  filter(!is.na(P) & P > 0) %>%
  mutate(
    CHR = as.integer(CHR),
    POS = as.numeric(POS),
    logP = -log10(P)
  )

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

set.seed(42)
gwas_sig <- gwas_plot_data %>% filter(logP >= 5)
gwas_nonsig <- gwas_plot_data %>% filter(logP < 5) %>% sample_frac(0.05)
gwas_plot_sample <- bind_rows(gwas_sig, gwas_nonsig)

panel_a <- ggplot(gwas_plot_sample, aes(x = cum_pos, y = logP, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_hline(yintercept = -log10(5e-8), color = "red", linetype = "dashed", linewidth = 0.5) +
  geom_hline(yintercept = -log10(1e-5), color = "cornflowerblue", linetype = "dashed", linewidth = 0.5) +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60")) +
  scale_x_continuous(
    breaks = chr_lengths$chr_mid,
    labels = chr_lengths$CHR,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title = "A",
    subtitle = "GWAS SNP-Level Association",
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5)) +
  coord_cartesian(ylim = c(0, max(gwas_plot_sample$logP) * 1.05))

magma_gene <- fread(combined_gene_file, header = TRUE)
if (!"P" %in% names(magma_gene) && "P_MULTI" %in% names(magma_gene)) {
  magma_gene$P <- magma_gene$P_MULTI
}
if (!all(c("CHR", "START", "STOP", "P") %in% names(magma_gene))) {
  stop("Combined MAGMA file is missing one of CHR/START/STOP/P.", call. = FALSE)
}

magma_gene <- magma_gene %>%
  mutate(
    CHR = as.integer(CHR),
    gene_mid = (START + STOP) / 2,
    logP = -log10(P)
  )

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

magma_fdr <- p.adjust(magma_gene$P, method = "BH")
magma_fdr_genes <- magma_gene[magma_fdr < 0.05, , drop = FALSE]
magma_fdr_threshold <- if (nrow(magma_fdr_genes)) -log10(max(magma_fdr_genes$P)) else NA_real_

panel_b <- ggplot(magma_plot_data, aes(x = cum_pos, y = logP, color = factor(CHR %% 2))) +
  geom_point(alpha = 0.6, size = 1) +
  geom_hline(yintercept = -log10(0.05 / nrow(magma_gene)), color = "red", linetype = "dashed", linewidth = 0.5) +
  {if (is.finite(magma_fdr_threshold)) geom_hline(yintercept = magma_fdr_threshold, color = "cornflowerblue", linetype = "dashed", linewidth = 0.5)} +
  scale_color_manual(values = c("0" = "gray30", "1" = "gray60")) +
  scale_x_continuous(
    breaks = chr_lengths_magma$chr_mid,
    labels = chr_lengths_magma$CHR,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(expand = c(0.02, 0)) +
  labs(
    title = "B",
    subtitle = "MAGMA Gene-Level Association",
    x = "Chromosome",
    y = expression(-log[10](p))
  ) +
  plot_theme +
  theme(plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5)) +
  coord_cartesian(ylim = c(0, max(magma_plot_data$logP) * 1.05))

pathway_results <- read.csv(catfish_csv, stringsAsFactors = FALSE)
final_p_col <- pick_final_p_col(pathway_results)
if (is.na(final_p_col)) {
  stop("Could not find a final omnibus p-value column in CATFISH CSV.", call. = FALSE)
}

top_pathways <- pathway_results %>%
  arrange(.data[[final_p_col]]) %>%
  slice_head(n = TOP_N_PATHWAYS) %>%
  mutate(
    pathway_label = paste0(pathway_name, " [", n_genes, "]"),
    pathway_label = factor(pathway_label, levels = rev(pathway_label))
  )

component_data <- top_pathways %>%
  select(
    pathway_label,
    ACAT = acat_p,
    Fisher = fisher_p,
    Stouffer = stouffer_p_analytic,
    minP = minp_p_analytic,
    TFisher = tfisher_p_analytic,
    Omnibus = all_of(final_p_col)
  ) %>%
  tidyr::pivot_longer(
    cols = c(ACAT, Fisher, Stouffer, minP, TFisher, Omnibus),
    names_to = "Method",
    values_to = "p_value"
  ) %>%
  mutate(
    logP = -log10(pmax(p_value, 1e-16)),
    Method = factor(Method, levels = c("ACAT", "Fisher", "Stouffer", "minP", "TFisher", "Omnibus"))
  )

panel_c <- ggplot(component_data, aes(x = Method, y = pathway_label, fill = logP)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "gray95",
    mid = "steelblue",
    high = "darkred",
    midpoint = 5,
    limits = c(0, max(12, ceiling(max(component_data$logP, na.rm = TRUE)))),
    name = expression(-log[10](p))
  ) +
  labs(
    title = "C",
    subtitle = "CATFISH Pathway Enrichment",
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
    plot.margin = margin(10, 10, 10, 10),
    plot.subtitle = element_text(size = 15, face = "bold", hjust = 0.5)
  )

fig <- (panel_a | panel_b) / panel_c + plot_layout(heights = c(1, 1.45))

png_file <- file.path(RESULTS_DIR, "Figure_SAP_Stem_volume_GWAS_MAGMA_CATFISH.png")
pdf_file <- file.path(RESULTS_DIR, "Figure_SAP_Stem_volume_GWAS_MAGMA_CATFISH.pdf")

ggsave(png_file, fig, width = 13.5, height = 12, dpi = 300, bg = "white")
ggsave(pdf_file, fig, width = 13.5, height = 12, bg = "white")

cat("Saved figure:\n")
cat("  ", png_file, "\n", sep = "")
cat("  ", pdf_file, "\n", sep = "")
