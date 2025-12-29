## ============================================================
## Arabidopsis (single group) version of your Fly figures:
##   1) Top-20 bubble plot (no sex facet)
##   2) UpSet overlap across 5 tests + Omni (top-N or p-cutoff)
##   3) Omni vs each method (2-circle Venn panels) within Arabidopsis
## ============================================================

library(tidyverse)
library(forcats)
library(viridis)
library(ComplexUpset)
library(patchwork)
library(ggplot2)

# -----------------------
# INPUT
# -----------------------
at_cold <- "magcat_omnibus_results_Arabidopsis/omni_minp_mvn.csv"
at_raw  <- readr::read_csv(at_cold, show_col_types = FALSE)

# -----------------------
# Helpers
# -----------------------
safe_p <- function(p, min_p = 1e-300) {
  p <- suppressWarnings(as.numeric(p))
  p[is.na(p)] <- 1
  p <- pmin(pmax(p, min_p), 1)
  p
}

pick_first_present <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (length(hit) == 0) NA_character_ else hit[1]
}

# Prefer omni_p_final if present; otherwise fall back
omni_col <- pick_first_present(at_raw, c("omni_p_final", "omni_p_mvn", "omni_p"))
if (is.na(omni_col)) stop("No omnibus p-value column found (expected omni_p_final / omni_p_mvn / omni_p).")

# Method columns (with small fallbacks)
method_candidates <- list(
  ACAT             = c("acat_p"),
  Fisher           = c("fisher_p"),
  Adaptive_TFisher = c("tfisher_p_analytic", "tfisher_p"),
  minP             = c("minp_p_analytic", "minp_p"),
  Stouffer         = c("stouffer_p_analytic", "stouffer_p")
)

method_cols <- purrr::map_chr(method_candidates, ~ pick_first_present(at_raw, .x))
method_cols <- method_cols[!is.na(method_cols)]  # drop missing methods safely

# -----------------------
# Theme (your style)
# -----------------------
plot_theme <- theme_minimal(base_size = 24) +
  theme(
    plot.title     = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
    axis.title.x   = element_text(size = 24, face = "bold"),
    axis.title.y   = element_text(size = 24, face = "bold"),
    axis.text.x    = element_text(size = 24, color = "black"),
    axis.text.y    = element_text(size = 24, color = "black"),
    axis.line      = element_line(color = "black"),
    panel.grid.major.x = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    legend.position = "top",
    legend.title    = element_blank(),
    legend.text     = element_text(size = 16),
    plot.margin     = margin(15, 15, 15, 15)
  )

# ============================================================
# 1) Top-20 bubble plot (single panel)
# ============================================================
top_n <- 20

at_top <- at_raw %>%
  mutate(
    omni_p_used = safe_p(.data[[omni_col]]),
    mlog10p     = -log10(omni_p_used),
    n_genes     = suppressWarnings(as.numeric(n_genes)),
    pathway_name = as.character(pathway_name)
  ) %>%
  arrange(omni_p_used) %>%
  slice_head(n = top_n) %>%
  mutate(pathway_name_ord = fct_reorder(pathway_name, mlog10p, .desc = FALSE))

p_at_bubble <- ggplot(at_top, aes(x = mlog10p, y = pathway_name_ord)) +
  geom_point(aes(size = n_genes, color = mlog10p), alpha = 0.95) +
  scale_color_viridis_c(option = "viridis", end = 0.98) +
  guides(
    size  = guide_legend(override.aes = list(alpha = 1)),
    color = guide_colorbar(barwidth = 14, barheight = 0.8)
  ) +
  labs(
    x     = expression(-log[10]("omni_p")),
    y     = NULL,
    size  = "n_genes",
    color = expression(-log[10]("omni_p")),
    title = paste0("Arabidopsis: Top ", top_n, " pathways (size = n_genes; color = -log10(omni))")
  ) +
  theme_bw() +
  plot_theme

quartz()
print(p_at_bubble)
ggsave("Figures/Fig_top15_bubbles_Arabidopsis.png", p_at_bubble,
       width = 18, height = 10, dpi = 300, bg = "white")


# ============================================================
# 2) UpSet: overlap across methods + Omni (single group)
#    Choose ONE mode:
#       A) top-N per method (default)
#       B) p-cutoff per method (optional)
# ============================================================

# ---------- A) top-N mode ----------
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

# Build mapping: column -> label
tests6 <- c(method_cols, setNames(omni_col, "Omni_MVN"))
tests6_named <- setNames(
  names(tests6),  # labels
  unname(tests6)  # columns
)
# Need names(cols_named)=columns, values=labels
cols_named <- setNames(names(tests6_named), names(tests6_named))
names(cols_named) <- names(tests6_named) # labels temporarily
# Fix properly:
cols_named <- setNames(names(tests6), unname(tests6))  # WRONG
# Let's do it clean:
cols_named <- setNames(
  object = c(names(method_cols), "Omni_MVN"),  # labels
  nm     = c(unname(method_cols), omni_col)    # columns
)

at_sets6 <- top_ids_by_col(at_raw, cols_named, top_n = 15, id_col = "pathway_id")
at_mem6  <- sets_to_df(at_sets6, id_col = "pathway_id")

p_up_at <- ComplexUpset::upset(
  at_mem6,
  intersect = names(at_sets6),
  name = "Arabidopsis top15 overlap (5 tests + Omni)",
  base_annotations = list("Intersection size" = intersection_size()),
  set_sizes = FALSE,
  themes = upset_default_themes(text = element_text(size = 26, color = "black"))
) + plot_theme +
  theme(axis.text.x = element_blank())

quartz()
print(p_up_at)
ggsave("Figures/Fig_Top15_overlap_Arabidopsis.png", p_up_at,
       width = 12, height = 8, dpi = 300, bg = "white")


# ---------- B) p-cutoff mode (optional) ----------
# p_cutoff <- 0.05
# use_fdr  <- FALSE
#
# ids_by_cutoff <- function(df, cols_named, p_cutoff = 0.05, id_col = "pathway_id", use_fdr = FALSE) {
#   out <- list()
#   for (col in names(cols_named)) {
#     lab <- cols_named[[col]]
#     tmp <- df %>% mutate(.p = safe_p(.data[[col]]))
#     if (use_fdr) tmp <- tmp %>% mutate(.p = p.adjust(.p, method = "BH"))
#     out[[lab]] <- tmp %>%
#       filter(.p <= p_cutoff) %>%
#       pull(.data[[id_col]]) %>%
#       unique()
#   }
#   out
# }
#
# at_sets6_cut <- ids_by_cutoff(at_raw, cols_named, p_cutoff = p_cutoff, use_fdr = use_fdr)
# at_mem6_cut  <- sets_to_df(at_sets6_cut, id_col = "pathway_id")
#
# p_up_at_cut <- ComplexUpset::upset(
#   at_mem6_cut,
#   intersect = names(at_sets6_cut),
#   name = paste0("Arabidopsis overlap (", ifelse(use_fdr, "BH q", "p"), " <= ", p_cutoff, "; 5 tests + Omni)"),
#   base_annotations = list("Intersection size" = intersection_size()),
#   set_sizes = FALSE,
#   themes = upset_default_themes(text = element_text(size = 26, color = "black"))
# ) + plot_theme + theme(axis.text.x = element_blank())
#
# print(p_up_at_cut)


# ============================================================
# 3) Omni vs each method (2-circle Venn panels) within Arabidopsis
#    (same geometry you used; single row of 5 panels)
# ============================================================

get_top_ids <- function(df, p_col, top_n = 100, id_col = "pathway_id") {
  df %>%
    mutate(.p = safe_p(.data[[p_col]])) %>%
    arrange(.p) %>%
    slice_head(n = top_n) %>%
    pull(.data[[id_col]]) %>%
    unique()
}

venn_counts <- function(method_ids, omni_ids) {
  n12 <- length(intersect(method_ids, omni_ids))
  n1  <- length(setdiff(method_ids, omni_ids))
  n2  <- length(setdiff(omni_ids, method_ids))
  list(n1 = n1, n12 = n12, n2 = n2)
}

circle_pts <- function(cx, cy, r = 1, n = 720) {
  t <- seq(0, 2*pi, length.out = n)
  tibble(x = cx + r*cos(t), y = cy + r*sin(t))
}
inside_circle <- function(x, y, cx, cy, r = 1) (x - cx)^2 + (y - cy)^2 <= r^2 + 1e-12
lens_polygon <- function(c1, c2, r = 1, n = 720) {
  p1 <- circle_pts(c1[1], c1[2], r, n) %>% filter(inside_circle(x, y, c2[1], c2[2], r))
  p2 <- circle_pts(c2[1], c2[2], r, n) %>% filter(inside_circle(x, y, c1[1], c1[2], r))
  P <- bind_rows(p1, p2)
  if (nrow(P) < 3) return(NULL)
  cx <- mean(P$x); cy <- mean(P$y)
  P %>% mutate(a = atan2(y - cy, x - cx)) %>% arrange(a) %>% select(x, y)
}

venn_two_custom <- function(n_method_only, n_overlap, n_omni_only,
                            method_label,
                            base_col = "#440154FF", overlap_col = "#FDE725FF",
                            r = 1, sep = 1.15) {

  c_method <- c(0,  sep/2)
  c_omni   <- c(0, -sep/2)

  top_circle <- circle_pts(c_method[1], c_method[2], r)
  bot_circle <- circle_pts(c_omni[1],   c_omni[2],   r)
  lens <- lens_polygon(c_method, c_omni, r = r)

  gg <- ggplot() +
    geom_polygon(data = top_circle, aes(x, y), fill = base_col, color = NA) +
    geom_polygon(data = bot_circle, aes(x, y), fill = base_col, color = NA)

  if (!is.null(lens)) gg <- gg + geom_polygon(data = lens, aes(x, y), fill = overlap_col, color = NA)

  gg +
    geom_path(data = top_circle, aes(x, y), color = "black", linewidth = 1) +
    geom_path(data = bot_circle, aes(x, y), color = "black", linewidth = 1) +
    coord_equal() +
    theme_void(base_size = 22) +
    theme(plot.margin = margin(10, 10, 10, 10)) +
    annotate("text", x = 0, y =  (sep/2 + r + 0.25), label = method_label,
             fontface = "bold", size = 7) +
    annotate("text", x = 0, y = -(sep/2 + r + 0.25), label = "OMNI",
             fontface = "bold", size = 7) +
    annotate("text", x = 0, y =  (sep/2 + 0.35), label = n_method_only, size = 7) +
    annotate("text", x = 0, y =  0,             label = n_overlap,     size = 7) +
    annotate("text", x = 0, y = -(sep/2 + 0.35), label = n_omni_only,   size = 7)
}

# top-N used for set definitions in these Venns
venn_top_n <- 15

omni_ids <- get_top_ids(at_raw, omni_col, top_n = venn_top_n)

p_venn_row <- purrr::imap(method_cols, function(col, lab) {
  method_ids <- get_top_ids(at_raw, col, top_n = venn_top_n)
  cc <- venn_counts(method_ids, omni_ids)
  venn_two_custom(cc$n1, cc$n12, cc$n2,
                  method_label = lab,
                  base_col = "#440154FF",      # purple circles
                  overlap_col = "#FDE725FF")   # yellow overlap
}) %>%
  purrr::reduce(`|`)

quartz()
print(p_venn_row)
ggsave("Figures/Fig_Arabidopsis_omni_vs_methods_venns.png", p_venn_row,
       width = 20, height = 6, dpi = 300, bg = "white")


# ============================================================
# Quick combined view (optional)
# ============================================================
# combined <- p_at_bubble / p_up_at / p_venn_row + plot_layout(heights = c(1.2, 1, 0.7))
# print(combined)
# ggsave("Figures/Fig_Arabidopsis_combo.png", combined, width = 18, height = 16, dpi = 300, bg = "white")
