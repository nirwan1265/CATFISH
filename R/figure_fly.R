library(tidyverse)
library(ggrepel)

f_female <- "magcat_omnibus_results_Fly/omni_minp_mvn_female.csv"
f_male   <- "magcat_omnibus_results_Fly/omni_minp_mvn_male.csv"

colnames(f_female)

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

p1 <- ggplot(top_bubbles, aes(x = mlog10p, y = pathway_name_ord)) +
  geom_point(aes(size = n_genes, color = mlog10p, shape = hit_type), alpha = 0.95) +
  facet_wrap(~sex, scales = "free_y") +
  scale_shape_manual(values = c("Female-only" = 17, "Male-only" = 15, "Shared" = 16)) +
  scale_color_viridis_c(option = "viridis", end = 0.98) +
  labs(
    x = expression(-log[10]("omni_p_final")),
    y = NULL,
    size = "n_genes",
    color = expression(-log[10]("omni_p_final")),
    shape = NULL,
    title = "Top 20 pathways per sex (size = n_genes; color = -log10(omni_p_final))"
  ) +
  theme_bw(base_size = 16) +
  theme(strip.text = element_text(face = "bold"))
  
quartz()
p1

# -----------------------
# Plot 2: Male vs Female scatter on union(top20)
#   - missing in one sex => p=1 => -log10(p)=0
#   - color = strongest signal in either sex (max of -log10 p)
#   - shape = hit_type, size = n_genes
# -----------------------
union20 <- bind_rows(
  topF %>% select(pathway_id, pathway_name, n_genes_F = n_genes, pF = omni_p_final, mF = mlog10p),
  topM %>% select(pathway_id, pathway_name, n_genes_M = n_genes, pM = omni_p_final, mM = mlog10p)
) %>%
  group_by(pathway_id) %>%
  summarise(
    pathway_name = dplyr::first(na.omit(pathway_name)),
    n_genes_F = dplyr::first(n_genes_F),
    n_genes_M = dplyr::first(n_genes_M),
    pF = dplyr::first(pF),
    pM = dplyr::first(pM),
    .groups = "drop"
  ) %>%
  mutate(
    pF = safe_p(pF),
    pM = safe_p(pM),
    mF = -log10(pF),
    mM = -log10(pM),
    m_any = pmax(mF, mM),
    n_genes = coalesce(n_genes_F, n_genes_M),
    hit_type = case_when(
      pathway_id %in% overlap_ids ~ "Shared",
      pathway_id %in% topF$pathway_id ~ "Female-only",
      pathway_id %in% topM$pathway_id ~ "Male-only",
      TRUE ~ "Other"
    )
  )


p2 <- ggplot(union20, aes(x = mM, y = mF)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2) +
  geom_point(aes(size = n_genes, color = m_any, shape = hit_type), alpha = 0.95) +
  scale_shape_manual(values = c("Female-only" = 17, "Male-only" = 15, "Shared" = 16, "Other" = 4)) +
  scale_color_viridis_c(option = "viridis", end = 0.98) +
  ggrepel::geom_text_repel(
    data = union20 %>% filter(hit_type == "Shared"),
    aes(label = pathway_name),
    max.overlaps = 50,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  labs(
    x = expression("Male  " * -log[10]("omni_p_final")),
    y = expression("Female  " * -log[10]("omni_p_final")),
    size = "n_genes",
    color = expression("max(-log[10](p))"),
    shape = NULL,
    title = "Male vs Female (union of top20): color = strongest signal in either sex"
  ) +
  theme_bw()

quartz()
p2


# Optional save:
ggsave("Figures/Fig_top20_bubbles_by_sex.png", p1, width = 20, height = 8, dpi = 300)
# ggsave("male_vs_female_union_top20.png", p2, width = 9, height = 7, dpi = 300)






#### VENN DIAGRAM
library(tidyverse)

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

female_sets <- top_ids_by_test(female_raw, top_n = 20)
male_sets   <- top_ids_by_test(male_raw,   top_n = 20)

# =========================================================
# 1) UpSet plot (RECOMMENDED) for 5 sets
# =========================================================
# install.packages("ComplexUpset")  # if needed
library(ComplexUpset)

sets_to_df <- function(sets_list) {
  # union of all ids
  all_ids <- sort(unique(unlist(sets_list)))
  # binary membership columns
  mem <- map_dfc(sets_list, ~ all_ids %in% .x) %>%
    mutate(pathway_id = all_ids, .before = 1)
  mem
}

female_mem <- sets_to_df(female_sets)
male_mem   <- sets_to_df(male_sets)

p_up_female <- ComplexUpset::upset(
  female_mem,
  intersect = names(female_sets),
  name = "Female top20 overlap (5 tests)"
)

p_up_male <- ComplexUpset::upset(
  male_mem,
  intersect = names(male_sets),
  name = "Male top20 overlap (5 tests)"
)

quartz()
p_up_female
p_up_male



library(tidyverse)
library(ComplexUpset)

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

top_ids_by_col <- function(df, cols_named, top_n = 20, id_col = "pathway_id") {
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

# ---- Female ----
female_sets6 <- top_ids_by_col(female_raw, tests6, top_n = 20, id_col = "pathway_id")
female_mem6  <- sets_to_df(female_sets6, id_col = "pathway_id")

p_up_female <- ComplexUpset::upset(
  female_mem6,
  intersect = names(female_sets6),
  name = "Female top20 overlap (5 tests + Omni_MVN)",
  base_annotations = list("Intersection size" = intersection_size())
)

# ---- Male ----
male_sets6 <- top_ids_by_col(male_raw, tests6, top_n = 20, id_col = "pathway_id")
male_mem6  <- sets_to_df(male_sets6, id_col = "pathway_id")

p_up_male <- ComplexUpset::upset(
  male_mem6,
  intersect = names(male_sets6),
  name = "Male top20 overlap (5 tests + Omni_MVN)",
  base_annotations = list("Intersection size" = intersection_size())
)

quartz()
p_up_female

quartz()
p_up_male



# =========================================================
# 2) Venn diagram (works, but cramped with 5)
# =========================================================
# install.packages("ggVennDiagram")  # if needed
library(ggVennDiagram)

p_venn_female <- ggVennDiagram(female_sets, label = "count") +
  ggtitle("Female: overlap of top 20 pathways across 5 tests")

p_venn_male <- ggVennDiagram(male_sets, label = "count") +
  ggtitle("Male: overlap of top 20 pathways across 5 tests")

quartz()
p_venn_female
quartz()
p_venn_male
