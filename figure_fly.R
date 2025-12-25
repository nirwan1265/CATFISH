library(tidyverse)
library(ggrepel)
library(tidyverse)
library(ComplexUpset)



# Pathway comparisong
f_female <- "magcat_omnibus_results_Fly/omni_minp_mvn_female.csv"
f_male   <- "magcat_omnibus_results_Fly/omni_minp_mvn_male.csv"

safe_p <- function(p, min_p = 1e-300) {
  p <- suppressWarnings(as.numeric(p))
  p[is.na(p)] <- 1
  p <- pmin(pmax(p, min_p), 1)
  p
}

prep_top <- function(df, sex, top_n = 20) {
  df %>%
    mutate(
      sex = sex,
      omni_p_final = safe_p(omni_p_final),
      mlog10p = -log10(omni_p_final),
      n_genes = suppressWarnings(as.numeric(n_genes))
    ) %>%
    arrange(omni_p_final) %>%
    slice_head(n = top_n) %>%
    select(sex, pathway_id, pathway_name, n_genes, omni_p_final, mlog10p)
}

female_raw <- readr::read_csv(f_female, show_col_types = FALSE)
male_raw   <- readr::read_csv(f_male,   show_col_types = FALSE)
colnames(female_raw)

topF <- prep_top(female_raw, "Female", top_n = 20)
topM <- prep_top(male_raw,   "Male",   top_n = 20)

overlap_ids <- intersect(topF$pathway_id, topM$pathway_id)

topF <- topF %>% mutate(hit_type = if_else(pathway_id %in% overlap_ids, "Shared", "Female-only"))
topM <- topM %>% mutate(hit_type = if_else(pathway_id %in% overlap_ids, "Shared", "Male-only"))

# -----------------------
# Plot 1: Faceted bubble dotplot (TOP 20 per sex)
#   color = -log10(p_final), size = n_genes, shape = hit_type
# -----------------------
top_bubbles <- bind_rows(topF, topM) %>%
  group_by(sex) %>%
  mutate(pathway_name_ord = fct_reorder(pathway_name, mlog10p, .desc = FALSE)) %>%
  ungroup()

plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_text(
      size = 24,      # X‐axis tick label size
      color = "black"
    ),
    axis.text.y    = element_text(
      size = 24,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),

    # ---- grids: light grey major grids for x + y, no minor grids ----
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),

    plot.margin    = margin(15, 15, 15, 15)

  )


p1 <- ggplot(top_bubbles, aes(x = mlog10p, y = pathway_name_ord)) +
  geom_point(aes(size = n_genes, color = mlog10p, shape = hit_type), alpha = 0.95) +
  facet_wrap(~sex, scales = "free_y") +
  scale_shape_manual(values = c("Female-only" = 17, "Male-only" = 15, "Shared" = 16)) +
  scale_color_viridis_c(option = "viridis", end = 0.98) +
  guides(shape = guide_legend(override.aes = list(size = 6)),
         size  = guide_legend(override.aes = list(alpha = 1)),
         color = guide_colorbar(barwidth = 14, barheight = 0.8)) +
  labs(
    x = expression(-log[10]("omni_p_final")),
    y = NULL,
    size = "n_genes",
    color = expression(-log[10]("omni_p_final")),
    shape = NULL,
    title = "Top 20 pathways per sex (size = n_genes; color = -log10(omni_p_final))"
  ) +
  theme_bw() +
  theme(strip.text = element_text(face = "bold")) +
  plot_theme

#quartz()
p1


# -----------------------
# Plot 2: Male vs Female scatter on union(top20)
#   - missing in one sex => p=1 => -log10(p)=0
#   - color = strongest signal in either sex (max of -log10 p)
#   - shape = hit_type, size = n_genes
# -----------------------
# union20 <- bind_rows(
#   topF %>% select(pathway_id, pathway_name, n_genes_F = n_genes, pF = omni_p_final, mF = mlog10p),
#   topM %>% select(pathway_id, pathway_name, n_genes_M = n_genes, pM = omni_p_final, mM = mlog10p)
# ) %>%
#   group_by(pathway_id) %>%
#   summarise(
#     pathway_name = dplyr::first(na.omit(pathway_name)),
#     n_genes_F = dplyr::first(n_genes_F),
#     n_genes_M = dplyr::first(n_genes_M),
#     pF = dplyr::first(pF),
#     pM = dplyr::first(pM),
#     .groups = "drop"
#   ) %>%
#   mutate(
#     pF = safe_p(pF),
#     pM = safe_p(pM),
#     mF = -log10(pF),
#     mM = -log10(pM),
#     m_any = pmax(mF, mM),
#     n_genes = coalesce(n_genes_F, n_genes_M),
#     hit_type = case_when(
#       pathway_id %in% overlap_ids ~ "Shared",
#       pathway_id %in% topF$pathway_id ~ "Female-only",
#       pathway_id %in% topM$pathway_id ~ "Male-only",
#       TRUE ~ "Other"
#     )
#   )
#
#
# p2 <- ggplot(union20, aes(x = mM, y = mF)) +
#   geom_abline(slope = 1, intercept = 0, linetype = 2) +
#   geom_point(aes(size = n_genes, color = m_any, shape = hit_type), alpha = 0.95) +
#   scale_shape_manual(values = c("Female-only" = 17, "Male-only" = 15, "Shared" = 16, "Other" = 4)) +
#   scale_color_viridis_c(option = "viridis", end = 0.98) +
#   ggrepel::geom_text_repel(
#     data = union20 %>% filter(hit_type == "Shared"),
#     aes(label = pathway_name),
#     max.overlaps = 50,
#     box.padding = 0.3,
#     point.padding = 0.2
#   ) +
#   labs(
#     x = expression("Male  " * -log[10]("omni_p_final")),
#     y = expression("Female  " * -log[10]("omni_p_final")),
#     size = "n_genes",
#     color = expression("max(-log[10](p))"),
#     shape = NULL,
#     title = "Male vs Female (union of top20): color = strongest signal in either sex"
#   ) +
#   theme_bw() +
#   plot_theme
#
# quartz()
# p2


# Optional save:
#ggsave("Figures/Fig_top20_bubbles_by_sex.png", p1, width = 24, height = 12, dpi = 300, bg = "white")
# ggsave("male_vs_female_union_top20.png", p2, width = 9, height = 7, dpi = 300)




#### VENN DIAGRAM
# library(tidyverse)
#
tests <- c(
  acat_p             = "ACAT",
  fisher_p           = "Fisher",
  tfisher_p_analytic  = "Adaptive_TFisher",
  minp_p_analytic     = "minP",
  stouffer_p_analytic = "Stouffer"
)

safe_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 1
  pmin(pmax(x, 1e-300), 1)
}

top_ids_by_test <- function(df, top_n = 20, id_col = "pathway_id") {
  out <- list()
  for (col in names(tests)) {
    tmp <- df %>%
      mutate(.p = safe_p(.data[[col]])) %>%
      arrange(.p) %>%
      slice_head(n = top_n) %>%
      pull(.data[[id_col]]) %>%
      unique()
    out[[tests[[col]]]] <- tmp
  }
  out
}

female_sets <- top_ids_by_test(female_raw, top_n = 100)
male_sets   <- top_ids_by_test(male_raw,   top_n = 100)

# =========================================================
# 1) UpSet plot (RECOMMENDED) for 5 sets
# =========================================================


safe_p <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 1
  pmin(pmax(x, 1e-300), 1)
}

# 5 tests + add omni_p_mvn as the 6th "set"
tests6 <- c(
  acat_p             = "ACAT",
  fisher_p           = "Fisher",
  tfisher_p_analytic  = "Adaptive_TFisher",
  minp_p_analytic     = "minP",
  stouffer_p_analytic = "Stouffer",
  omni_p_mvn          = "Omni_MVN"
)

top_ids_by_col <- function(df, cols_named, top_n = 15, id_col = "pathway_id") {
  out <- list()
  for (col in names(cols_named)) {
    lab <- cols_named[[col]]
    ids <- df %>%
      mutate(.p = safe_p(.data[[col]])) %>%
      arrange(.p) %>%
      slice_head(n = top_n) %>%
      pull(.data[[id_col]]) %>%
      unique()
    out[[lab]] <- ids
  }
  out
}

sets_to_df <- function(sets_list, id_col = "pathway_id") {
  all_ids <- sort(unique(unlist(sets_list)))
  mem <- purrr::map_dfc(sets_list, ~ all_ids %in% .x)
  mem <- dplyr::mutate(mem, !!id_col := all_ids, .before = 1)
  mem
}

upset_big_grid <- theme_minimal(base_size = 26) +
  theme(
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor   = element_blank(),
    axis.text          = element_text(color = "black"),
    axis.title         = element_text(color = "black", face = "bold")
  )

# ---- Female ----
female_sets6 <- top_ids_by_col(female_raw, tests6, top_n = 10, id_col = "pathway_id")
female_mem6  <- sets_to_df(female_sets6, id_col = "pathway_id")


plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(
      size   = 14,
      face   = "bold",
      hjust  = 0.5,
      margin = margin(b = 10)
    ),
    axis.title.x   = element_text(
      size = 24,      # X‐axis title size
      face = "bold"
    ),
    axis.title.y   = element_text(
      size = 24,      # Y‐axis title size
      face = "bold"
    ),
    axis.text.x    = element_blank(),
    axis.text.y    = element_text(
      size = 24,      # Y‐axis tick label size
      color = "black"
    ),
    axis.line      = element_line(color = "black"),

    # ---- grids: light grey major grids for x + y, no minor grids ----
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),

    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),

    plot.margin    = margin(15, 15, 15, 15)

  )


p_up_female <- ComplexUpset::upset(
  female_mem6,
  intersect = names(female_sets6),
  name = "Female top10 overlap (5 tests + Omni_MVN)",
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes = FALSE,
  themes = upset_default_themes(
    text = element_text(size = 26, color = "black")  # <-- BIGGER TEXT,
  )
) + plot_theme



# ---- Male ----
male_sets6 <- top_ids_by_col(male_raw, tests6, top_n = 10, id_col = "pathway_id")
male_mem6  <- sets_to_df(male_sets6, id_col = "pathway_id")

p_up_male <- ComplexUpset::upset(
  male_mem6,
  intersect = names(male_sets6),
  name = "Male top10 overlap (5 tests + Omni_MVN)",
  base_annotations = list(
    "Intersection size" = intersection_size()
  ),
  set_sizes = FALSE,
  themes = upset_default_themes(
    text = element_text(size = 26, color = "black")  # <-- BIGGER TEXT
  )
) + plot_theme


quartz()
p_up_female

quartz()
p_up_male



# Save ggsave
ggsave("Figures/Fig_Top15_overlap_female.png", p_up_female, width = 12, height = 8, dpi = 300, bg = "white")
ggsave("Figures/Fig_Top15_overlap_male.png", p_up_male, width = 12, height = 8, dpi = 300, bg = "white")






# # =========================================================
# # 0) CONFIG
# #    - Replace top-N with p-value (or q-value) cutoff
# #    - Keeps SAME output structure for ComplexUpset
# # =========================================================
#
# # 5 tests + Omni_MVN as the 6th "set"
# tests6 <- c(
#   acat_p              = "ACAT",
#   fisher_p            = "Fisher",
#   tfisher_p_analytic   = "Adaptive_TFisher",
#   minp_p_analytic      = "minP",
#   stouffer_p_analytic  = "Stouffer",
#   omni_p_mvn           = "Omni_MVN"
# )
#
# safe_p <- function(x) {
#   x <- suppressWarnings(as.numeric(x))
#   x[is.na(x)] <- 1
#   pmin(pmax(x, 1e-300), 1)
# }
#
# # =========================================================
# # 1) BUILD SETS BY P-VALUE CUTOFF  (NO top-N)
# # =========================================================
#
# # Choose ONE:
# p_cutoff <- 0.1         # e.g., 0.001
# # p_cutoff <- 1e-2        # e.g., 0.01
# # p_cutoff <- 5e-2        # e.g., 0.05
#
# # If you want FDR instead of raw p, set use_fdr <- TRUE
# use_fdr <- FALSE         # TRUE = BH-adjust within each sex & column
#
# ids_by_cutoff <- function(df, cols_named, p_cutoff = 1e-3, id_col = "pathway_id", use_fdr = FALSE) {
#   out <- list()
#
#   for (col in names(cols_named)) {
#     lab <- cols_named[[col]]
#
#     tmp <- df %>%
#       mutate(.p = safe_p(.data[[col]]))
#
#     if (use_fdr) {
#       tmp <- tmp %>% mutate(.p = p.adjust(.p, method = "BH"))
#     }
#
#     ids <- tmp %>%
#       filter(.p <= p_cutoff) %>%
#       pull(.data[[id_col]]) %>%
#       unique()
#
#     out[[lab]] <- ids
#   }
#
#   out
# }
#
# sets_to_df <- function(sets_list, id_col = "pathway_id") {
#   all_ids <- sort(unique(unlist(sets_list)))
#   mem <- purrr::map_dfc(sets_list, ~ all_ids %in% .x)
#   mem <- dplyr::mutate(mem, !!id_col := all_ids, .before = 1)
#   mem
# }
#
# # =========================================================
# # 2) THEMES (UNCHANGED)
# # =========================================================
#
# upset_big_grid <- theme_minimal(base_size = 26) +
#   theme(
#     panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
#     panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
#     panel.grid.minor   = element_blank(),
#     axis.text          = element_text(color = "black"),
#     axis.title         = element_text(color = "black", face = "bold")
#   )
#
# plot_theme <- theme_minimal(base_size = 24) +
#   theme(
#     plot.title     = element_text(
#       size   = 14,
#       face   = "bold",
#       hjust  = 0.5,
#       margin = margin(b = 10)
#     ),
#     axis.title.x   = element_text(
#       size = 24,
#       face = "bold"
#     ),
#     axis.title.y   = element_text(
#       size = 24,
#       face = "bold"
#     ),
#     axis.text.x    = element_blank(),
#     axis.text.y    = element_text(
#       size = 24,
#       color = "black"
#     ),
#     axis.line      = element_line(color = "black"),
#
#     panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
#     panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
#     panel.grid.minor.x = element_blank(),
#     panel.grid.minor.y = element_blank(),
#
#     legend.position = "top",
#     legend.title    = element_blank(),
#     legend.text     = element_text(size = 16),
#
#     plot.margin    = margin(15, 15, 15, 15)
#   )
#
# # =========================================================
# # 3) FEMALE UPSET (P-CUTOFF)
# # =========================================================
#
# female_sets6 <- ids_by_cutoff(
#   female_raw,
#   tests6,
#   p_cutoff = p_cutoff,
#   id_col   = "pathway_id",
#   use_fdr  = use_fdr
# )
#
# female_mem6 <- sets_to_df(female_sets6, id_col = "pathway_id")
#
# female_title <- if (use_fdr) {
#   paste0("Female overlap (BH q <= ", format(p_cutoff, scientific = TRUE), "; 5 tests + Omni_MVN)")
# } else {
#   paste0("Female overlap (p <= ", format(p_cutoff, scientific = TRUE), "; 5 tests + Omni_MVN)")
# }
#
# p_up_female <- ComplexUpset::upset(
#   female_mem6,
#   intersect = names(female_sets6),
#   name = female_title,
#   base_annotations = list(
#     "Intersection size" = ComplexUpset::intersection_size()
#   ),
#   set_sizes = FALSE,
#   themes = ComplexUpset::upset_default_themes(
#     text = element_text(size = 26, color = "black")
#   )
# ) + plot_theme
#
# # =========================================================
# # 4) MALE UPSET (P-CUTOFF)
# # =========================================================
#
# male_sets6 <- ids_by_cutoff(
#   male_raw,
#   tests6,
#   p_cutoff = p_cutoff,
#   id_col   = "pathway_id",
#   use_fdr  = use_fdr
# )
#
# male_mem6 <- sets_to_df(male_sets6, id_col = "pathway_id")
#
# male_title <- if (use_fdr) {
#   paste0("Male overlap (BH q <= ", format(p_cutoff, scientific = TRUE), "; 5 tests + Omni_MVN)")
# } else {
#   paste0("Male overlap (p <= ", format(p_cutoff, scientific = TRUE), "; 5 tests + Omni_MVN)")
# }
#
# p_up_male <- ComplexUpset::upset(
#   male_mem6,
#   intersect = names(male_sets6),
#   name = male_title,
#   base_annotations = list(
#     "Intersection size" = ComplexUpset::intersection_size()
#   ),
#   set_sizes = FALSE,
#   themes = ComplexUpset::upset_default_themes(
#     text = element_text(size = 26, color = "black")
#   )
# ) + plot_theme
#
# # =========================================================
# # 5) VIEW + SAVE (UNCHANGED)
# # =========================================================
#
# quartz()
# p_up_female
#
# quartz()
# p_up_male
#
# # Use cutoff in filenames so you don't overwrite figures
# suffix <- if (use_fdr) {
#   paste0("q", gsub("\\+|\\-", "", format(p_cutoff, scientific = TRUE)))
# } else {
#   paste0("p", gsub("\\+|\\-", "", format(p_cutoff, scientific = TRUE)))
# }
#
# ggsave(paste0("Figures/Fig_overlap_female_", suffix, ".png"),
#        p_up_female, width = 12, height = 8, dpi = 300, bg = "white")
#
# ggsave(paste0("Figures/Fig_overlap_male_", suffix, ".png"),
#        p_up_male, width = 12, height = 8, dpi = 300, bg = "white")




# =========================================================
# 2) Venn diagram (works, but cramped with 5)
# =========================================================
# install.packages("ggVennDiagram")  # if needed
# library(ggVennDiagram)
#
# p_venn_female <- ggVennDiagram(female_sets, label = "count") +
#   ggtitle("Female: overlap of top 20 pathways across 5 tests")
#
# p_venn_male <- ggVennDiagram(male_sets, label = "count") +
#   ggtitle("Male: overlap of top 20 pathways across 5 tests")
#
# quartz()
# p_venn_female
# quartz()
# p_venn_male

male_mem6
female_mem6

# ---- helper: extract pathway IDs where a given method column is TRUE ----
get_ids <- function(mem_tbl, method_col) {
  mem_tbl %>%
    filter(.data[[method_col]] %in% TRUE) %>%
    pull(pathway_id) %>%
    unique()
}

# ---- one Venn in *your* exact style ----
venn_one_style <- function(method_col, title_txt) {
  male_ids   <- get_ids(male_mem6, method_col)
  female_ids <- get_ids(female_mem6, method_col)

  venn_list <- list(Male = male_ids, Female = female_ids)

  ggVennDiagram(
    venn_list,
    label       = "count",
    label_alpha = 0,
    label_size  = 8
  ) +
    scale_fill_viridis_c(option = "viridis", direction = 1) +
    theme_void(base_size = 24) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 22, face = "bold", hjust = 0.5)
    ) +
    ggtitle(title_txt)
}

# ---- 6 methods ----
p_acat     <- venn_one_style("ACAT",             "ACAT")
p_fisher   <- venn_one_style("Fisher",           "Fisher")
p_tfisher  <- venn_one_style("Adaptive_TFisher", "Adaptive TFisher")
p_minp     <- venn_one_style("minP",             "minP")
p_stouffer <- venn_one_style("Stouffer",         "Stouffer")
p_omni     <- venn_one_style("Omni_MVN",         "Omni MVN")

# ---- grid 3x2 ----
p_grid <- (p_acat | p_fisher | p_tfisher) /
  (p_minp | p_stouffer | p_omni)

quartz()
p_grid

# optional save:
ggsave("Venn_male_vs_female_pathways.png", p_grid, width = 12, height = 8, dpi = 300)




# =========================================================
### #INDIVIDUAL GWAS
# =========================================================
male_gwas=read.csv("Tables/Individual_GWAS_male_starvation.csv")
female_gwas=read.csv("Tables/Individual_GWAS_female_starvation.csv")

head(male_gwas)
head(female_gwas)

library(tidyverse)
library(ggVennDiagram)

# Assume you already have: male_gwas, female_gwas
# We’ll take TOP 50 by smallest P.value (unique genes)

topN_genes <- function(df, gene_col = "GeneID", p_col = "P.value", top_n = 50) {
  df %>%
    mutate(
      .p = suppressWarnings(as.numeric(.data[[p_col]])),
      .gene = as.character(.data[[gene_col]])
    ) %>%
    filter(!is.na(.p), !is.na(.gene), .gene != "") %>%
    # if gene IDs have a "gene:" prefix (female does), strip it so IDs match
    mutate(.gene = sub("^gene:", "", .gene)) %>%
    arrange(.p) %>%
    distinct(.gene, .keep_all = TRUE) %>%   # keep best SNP per gene
    slice_head(n = top_n) %>%
    pull(.gene) %>%
    unique()
}

male_top50   <- topN_genes(male_gwas,   gene_col = "GeneID", p_col = "P.value", top_n = 50)
female_top50 <- topN_genes(female_gwas, gene_col = "GeneID", p_col = "P.value", top_n = 50)

venn_list_top50 <- list(
  Male   = male_top50,
  Female = female_top50
)

p_venn_top50 <- ggVennDiagram(
  venn_list_top50,
  label = "count",
  label_alpha = 0,
  label_size = 8
) +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  theme_void(base_size = 24) +
  theme(legend.position = "none",
        plot.title     = element_blank()) +
  ggtitle("Top 50 GWAS genes overlap (by smallest P.value)") +


quartz()
p_venn_top50

# optional: counts + shared IDs
cat("Male top50:", length(male_top50), "\n")
cat("Female top50:", length(female_top50), "\n")
cat("Shared:", length(intersect(male_top50, female_top50)), "\n")

# Save the file
ggsave("Figures/Fig_Top50_Individual_GWAS_genes_venn.png", p_venn_top50, width = 8, height = 6, dpi = 300, bg = "white")





# =========================================================
### MAGMA results
# =========================================================
male_magma <- read.csv("Tables/MAGMA_genes_all_raw_male_fly.csv")
female_magma <- read.csv("Tables/MAGMA_genes_all_raw_female_fly.csv")
head(male_magma)
head(female_magma)


topN_magma_genes <- function(df, top_n = 50) {
  df %>%
    mutate(
      P = suppressWarnings(as.numeric(P)),
      GENE = as.character(GENE)
    ) %>%
    filter(!is.na(P), !is.na(GENE), GENE != "") %>%
    arrange(P) %>%
    distinct(GENE, .keep_all = TRUE) %>%   # one row per gene (best P)
    slice_head(n = top_n) %>%
    pull(GENE) %>%
    unique()
}

male_top50_magma   <- topN_magma_genes(male_magma,   top_n = 50)
female_top50_magma <- topN_magma_genes(female_magma, top_n = 50)

venn_list_top50_magma <- list(
  Male   = male_top50_magma,
  Female = female_top50_magma
)

p_venn_top50_magma <- ggVennDiagram(
  venn_list_top50_magma,
  label       = "count",
  label_alpha = 0,
  label_size  = 8
) +
  scale_fill_viridis_c(option = "viridis", direction = 1) +
  theme_void(base_size = 24) +
  theme(
    legend.position = "none",
    plot.title      = element_text(size = 22, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Top 50 MAGMA genes overlap (by smallest P)")

quartz()
p_venn_top50_magma

# optional counts
cat("Male top50:", length(male_top50_magma), "\n")
cat("Female top50:", length(female_top50_magma), "\n")
cat("Shared:", length(intersect(male_top50_magma, female_top50_magma)), "\n")

# Save the file
ggsave("Figures/Fig_Top50_genes_venn_MAGMA.png", p_venn_top50_magma, width = 8, height = 6, dpi = 300, bg = "white")



##### OMNI VS CATFISH
library(tidyverse)
library(patchwork)

# -------------------------
# Helpers: set extraction + counts
# -------------------------
get_ids <- function(mem_tbl, method_col) {
  mem_tbl %>%
    filter(.data[[method_col]] %in% TRUE) %>%
    pull(pathway_id) %>%
    unique()
}

venn_counts <- function(method_ids, omni_ids) {
  n12 <- length(intersect(method_ids, omni_ids))
  n1  <- length(setdiff(method_ids, omni_ids))  # method-only
  n2  <- length(setdiff(omni_ids, method_ids))  # omni-only
  list(n1 = n1, n12 = n12, n2 = n2)
}

# -------------------------
# Geometry: draw 2-circle Venn (vertical) with a TRUE overlap lens polygon
# -------------------------
circle_pts <- function(cx, cy, r = 1, n = 720) {
  t <- seq(0, 2*pi, length.out = n)
  tibble(x = cx + r*cos(t), y = cy + r*sin(t))
}

inside_circle <- function(x, y, cx, cy, r = 1) {
  (x - cx)^2 + (y - cy)^2 <= r^2 + 1e-12
}

lens_polygon <- function(c1, c2, r = 1, n = 720) {
  # points on circle1 that lie inside circle2 + points on circle2 that lie inside circle1
  p1 <- circle_pts(c1[1], c1[2], r, n) %>%
    filter(inside_circle(x, y, c2[1], c2[2], r))
  p2 <- circle_pts(c2[1], c2[2], r, n) %>%
    filter(inside_circle(x, y, c1[1], c1[2], r))

  P <- bind_rows(p1, p2)
  if (nrow(P) < 3) return(NULL)

  cx <- mean(P$x); cy <- mean(P$y)
  P <- P %>%
    mutate(a = atan2(y - cy, x - cx)) %>%
    arrange(a) %>%
    select(x, y)
  P
}

venn_two_custom <- function(n_method_only, n_overlap, n_omni_only,
                            method_label,
                            base_col, overlap_col,
                            r = 1, sep = 1.15) {

  # vertical centers: method (top), omni (bottom)
  c_method <- c(0,  sep/2)
  c_omni   <- c(0, -sep/2)

  top_circle <- circle_pts(c_method[1], c_method[2], r)
  bot_circle <- circle_pts(c_omni[1],   c_omni[2],   r)
  lens <- lens_polygon(c_method, c_omni, r = r)

  gg <- ggplot() +
    # fills
    geom_polygon(data = top_circle, aes(x, y), fill = base_col, color = NA) +
    geom_polygon(data = bot_circle, aes(x, y), fill = base_col, color = NA)

  if (!is.null(lens)) {
    gg <- gg + geom_polygon(data = lens, aes(x, y), fill = overlap_col, color = NA)
  }

  # outlines
  gg <- gg +
    geom_path(data = top_circle, aes(x, y), color = "black", linewidth = 1) +
    geom_path(data = bot_circle, aes(x, y), color = "black", linewidth = 1) +
    coord_equal() +
    theme_void(base_size = 22) +
    theme(plot.margin = margin(10, 10, 10, 10))

  # labels (TOP = method, BOTTOM = OMNI)
  gg <- gg +
    annotate("text", x = 0, y =  (sep/2 + r + 0.25), label = method_label,
             fontface = "bold", size = 7) +
    annotate("text", x = 0, y = -(sep/2 + r + 0.25), label = "OMNI",
             fontface = "bold", size = 7)

  # counts (top-only / overlap / bottom-only)
  gg <- gg +
    annotate("text", x = 0, y =  (sep/2 + 0.35), label = n_method_only, size = 7) +
    annotate("text", x = 0, y =  0,             label = n_overlap,     size = 7) +
    annotate("text", x = 0, y = -(sep/2 + 0.35), label = n_omni_only,   size = 7)

  gg
}

# -------------------------
# Colors you asked for
# Female: purple circles, yellow overlap
# Male:   yellow circles, purple overlap
# -------------------------
female_base    <- "#440154FF"  # purple
female_overlap <- "#FDE725FF"  # yellow
male_base      <- "#FDE725FF"  # yellow
male_overlap   <- "#440154FF"  # purple

# -------------------------
# Build 5 panels per sex: OMNI vs each method
# -------------------------
make_row <- function(mem_tbl, sex = c("female", "male")) {
  sex <- match.arg(sex)

  omni_ids <- get_ids(mem_tbl, "Omni_MVN")

  methods <- c(
    ACAT             = "ACAT",
    Fisher           = "Fisher",
    Adaptive_TFisher = "TFisher",
    minP             = "minP",
    Stouffer         = "Stouffer"
  )

  plots <- lapply(names(methods), function(col) {
    method_ids <- get_ids(mem_tbl, col)
    cc <- venn_counts(method_ids, omni_ids)

    if (sex == "female") {
      venn_two_custom(cc$n1, cc$n12, cc$n2, methods[[col]],
                      base_col = female_base, overlap_col = female_overlap)
    } else {
      venn_two_custom(cc$n1, cc$n12, cc$n2, methods[[col]],
                      base_col = male_base, overlap_col = male_overlap)
    }
  })

  plots[[1]] | plots[[2]] | plots[[3]] | plots[[4]] | plots[[5]]
}

grid_female <- make_row(female_mem6, "female")
grid_male   <- make_row(male_mem6,   "male")

total_grid <- grid_female / grid_male + plot_layout(heights = c(1, 1))

quartz()
total_grid

ggsave("venn_omni_vs_methods_female_male_fixedcolors.png",
       total_grid, width = 18, height = 8, dpi = 300)
