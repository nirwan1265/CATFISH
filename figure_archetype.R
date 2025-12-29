# Full script (I–VI) with Archetype VI (CEAB) INCLUDED.
# Also avoids slice() entirely (uses base indexing) to prevent Rle/list conflicts.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(scales)
})

set.seed(1)

# ---------------- helpers ----------------
rank_curve <- function(p, label, group = "in_set") {
  p <- p[is.finite(p) & !is.na(p)]
  p <- pmin(pmax(p, .Machine$double.xmin), 1)
  tibble(
    archetype = label,
    group     = group,
    rank      = seq_along(p) / length(p),
    p         = sort(p),
    mlog10p   = -log10(sort(p))
  )
}

theme_natureish <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title    = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size),
      axis.title    = element_text(size = base_size),
      axis.text     = element_text(size = base_size - 1, color = "black"),
      panel.border  = element_rect(fill = NA, linewidth = 0.35),
      plot.margin   = margin(8, 8, 8, 8),
      legend.position = "none"
    )
}

anno_box <- function(x, y, label, size = 3.6) {
  geom_label(
    inherit.aes = FALSE,
    data = data.frame(x = x, y = y, label = label),
    aes(x = x, y = y, label = label),
    hjust = 0, vjust = 1,
    size = size,
    label.size = 0,
    fill = scales::alpha("white", 0.85)
  )
}

add_thresholds <- function(g, x = 0.02, size = 3) {
  thr <- tibble(
    y   = -log10(c(5e-2, 1e-3, 2.5e-6)),
    lab = c("0.05", "1e-3", "2.5e-6")
  )
  g +
    geom_hline(data = thr, aes(yintercept = y), linewidth = 0.3, alpha = 0.55) +
    annotate("text", x = x, y = thr$y, label = thr$lab,
             hjust = 0, vjust = -0.3, size = size)
}

# ---------------- simulate archetype patterns ----------------
n <- 220
p_I   <- c(10^runif(6, -14, -7), runif(n - 6))                              # SDA
p_II  <- pmin(pmax(rbeta(n, 0.55, 2.6), 1e-8), 1)                           # CME
p_III <- pmin(pmax(rbeta(n, 0.95, 1.25), 1e-8), 1)                          # DPS
p_IV  <- c(10^runif(4, -12, -7), rbeta(40, 0.65, 2.8), runif(n - 44))       # HDS
p_V   <- c(10^runif(1, -20, -12), runif(n - 1))                             # SGP

# VI — CEAB: in-set shifted vs out-set background (competitive enrichment)
n_in  <- 120
n_out <- 400
p_VI_in  <- pmin(pmax(rbeta(n_in, 0.85, 1.9), 1e-8), 1)
p_VI_out <- runif(n_out)

# ---------------- build plotting data ----------------
df <- bind_rows(
  rank_curve(p_I,   "Archetype I — SDA"),
  rank_curve(p_II,  "Archetype II — CME"),
  rank_curve(p_III, "Archetype III — DPS"),
  rank_curve(p_IV,  "Archetype IV — HDS"),
  rank_curve(p_V,   "Archetype V — SGP")
)

df_vi <- bind_rows(
  rank_curve(p_VI_in,  "Archetype VI — CEAB", group = "in_set"),
  rank_curve(p_VI_out, "Archetype VI — CEAB", group = "out_set")
)

# ---------------- I — SDA ----------------
sda_dat <- df %>% filter(archetype == "Archetype I — SDA") %>% arrange(rank)
topM <- 6
sda_top <- sda_dat[seq_len(topM), , drop = FALSE]

k_1e3_sda <- sum(sda_dat$p < 1e-3)
k_gw_sda  <- sum(sda_dat$p < 2.5e-6)

p1 <- ggplot(sda_dat, aes(rank, mlog10p)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 0.55, alpha = 0.75) +
  geom_point(data = sda_top, size = 2) +
  coord_cartesian(ylim = c(0, 12.5)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype I — Sparse Driver (SDA)",
    subtitle = "several extreme genes; not just one",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p1 <- add_thresholds(p1) +
  anno_box(0.12, 12.35,
           paste0("highlighted top ", topM, " genes\n",
                  "K(p<1e-3)=", k_1e3_sda, "   K(p<2.5e-6)=", k_gw_sda))

# ---------------- II — CME ----------------
cme_dat <- df %>% filter(archetype == "Archetype II — CME") %>% arrange(rank)
p2 <- ggplot(cme_dat, aes(rank, mlog10p)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype II — Coordinated Moderate (CME)",
    subtitle = "broad moderate signal; no single extreme",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p2 <- add_thresholds(p2)

# ---------------- III — DPS ----------------
dps_dat <- df %>% filter(archetype == "Archetype III — DPS") %>% arrange(rank)
p3 <- ggplot(dps_dat, aes(rank, mlog10p)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype III — Diffuse Polygenic Shift (DPS)",
    subtitle = "slight upward shift; almost none cross strict thresholds",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p3 <- add_thresholds(p3)

# ---------------- IV — HDS ----------------
hds_dat <- df %>% filter(archetype == "Archetype IV — HDS") %>% arrange(rank)
p4 <- ggplot(hds_dat, aes(rank, mlog10p)) +
  geom_line(linewidth = 0.7) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype IV — Hybrid Driver–Support (HDS)",
    subtitle = "few drivers + a shoulder of moderate genes",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p4 <- add_thresholds(p4)

# ---------------- V — SGP ----------------
sgp_dat <- df %>% filter(archetype == "Archetype V — SGP") %>% arrange(rank)
sgp_tail <- sgp_dat[-1, , drop = FALSE]
sgp_tail$rank_tail    <- seq_along(sgp_tail$p) / length(sgp_tail$p)
sgp_tail$mlog10p_tail <- -log10(sgp_tail$p)

k_1e3_sgp <- sum(sgp_dat$p < 1e-3)
k_gw_sgp  <- sum(sgp_dat$p < 2.5e-6)
top1 <- sgp_dat[1, , drop = FALSE]

p5 <- ggplot() +
  geom_line(data = sgp_dat,  aes(rank, mlog10p), linewidth = 0.7) +
  geom_point(data = sgp_dat, aes(rank, mlog10p), size = 0.55, alpha = 0.75) +
  geom_line(data = sgp_tail, aes(rank_tail, mlog10p_tail), linewidth = 0.7, linetype = "22") +
  geom_point(data = top1, aes(rank, mlog10p), size = 2) +
  coord_cartesian(ylim = c(0, 8.5)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype V — Single-Gene Proxy (SGP)",
    subtitle = "one dominant gene; tail-only curve ~null after removing it",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p5 <- add_thresholds(p5) +
  anno_box(0.12, 8.35,
           paste0("K(p<1e-3)=", k_1e3_sgp, "   K(p<2.5e-6)=", k_gw_sgp, "\n",
                  "dashed: tail-only (top gene removed)"))

# ---------------- VI — CEAB (THIS is archetype 6) ----------------
p6 <- ggplot(df_vi, aes(rank, mlog10p, linetype = group)) +
  geom_line(linewidth = 0.75) +
  scale_linetype_manual(values = c(in_set = "solid", out_set = "22")) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype VI — Competitive Enrichment (CEAB)",
    subtitle = "in-set curve sits above background (βs > 0 in MAGMA competitive test)",
    x = "Gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()

p6 <- p6 +
  anno_box(0.08, max(df_vi$mlog10p) * 0.92,
           "solid: in-set\n dashed: out-of-set", size = 3.2)

# ---------------- compose 2x3 figure (includes VI) ----------------
fig <- (p1 + p2 + p3) / (p4 + p5 + p6) +
  plot_annotation(
    title = "Pathway signal archetypes as gene-level p-value rank curves",
    theme = theme(plot.title = element_text(face = "bold", size = 13))
  )

#quartz()
fig

# Save:
#ggsave("Fig_archetypesI.png", p1, width = 7.2, height = 4.6, units = "in", dpi = 600)
#ggsave("Fig_archetypesII.png", p2, width = 7.2, height = 4.6, units = "in", dpi = 600)
#ggsave("Fig_archetypesIII.png", p3, width = 7.2, height = 4.6, units = "in", dpi = 600)
#ggsave("Fig_archetypesIV.png", p4, width = 7.2, height = 4.6, units = "in", dpi = 600)
#ggsave("Fig_archetypesV.png", p5, width = 7.2, height = 4.6, units = "in", dpi = 600)
#ggsave("Fig_archetypesVI.png", p6, width = 7.2, height = 4.6, units = "in", dpi = 600)

# ggsave("Fig_archetypes_rankcurves.pdf", fig, width = 7.2, height = 4.6, units = "in", device = cairo_pdf)
# ggsave("Fig_archetypes_rankcurves.png", fig, width = 7.2, height = 4.6, units = "in", dpi = 600)
