# Archetype rank-curve figure for the five core pathway signal families:
# I SDA, II CME, III DPS, IV HDS, V SGP.
# Simulates gene-level Z from MVN: Z ~ N(mu, Sigma), then converts to p-values.

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
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

# ---------------- simulate archetype patterns from MVN ----------------
make_exchangeable_sigma <- function(m, rho) {
  matrix(rho, nrow = m, ncol = m) + diag(1 - rho, m)
}

simulate_mvn_p <- function(mu, sigma, one_sided = TRUE) {
  m <- length(mu)
  z <- as.numeric(mu + t(chol(sigma)) %*% rnorm(m))
  p <- if (one_sided) {
    pnorm(z, lower.tail = FALSE)
  } else {
    2 * pnorm(-abs(z))
  }
  pmin(pmax(p, .Machine$double.xmin), 1)
}

m <- 220
rho <- 0.20
sigma <- make_exchangeable_sigma(m = m, rho = rho)

# I — SDA: few strong drivers
mu_I <- rep(0, m)
idx_I <- sample.int(m, 6)
mu_I[idx_I] <- 3.2
p_I <- simulate_mvn_p(mu_I, sigma)

# II — CME: many moderate same-direction effects
mu_II <- rep(0, m)
idx_II <- sample.int(m, round(0.45 * m))
mu_II[idx_II] <- 1.35
p_II <- simulate_mvn_p(mu_II, sigma)

# III — DPS: diffuse weak shift across most genes
mu_III <- rep(0, m)
idx_III <- sample.int(m, round(0.85 * m))
mu_III[idx_III] <- 0.45
p_III <- simulate_mvn_p(mu_III, sigma)

# IV — HDS: few strong + broader moderate support
mu_IV <- rep(0, m)
idx_driver <- sample.int(m, 4)
idx_support <- sample(setdiff(seq_len(m), idx_driver), 40)
mu_IV[idx_driver] <- 3.0
mu_IV[idx_support] <- 1.1
p_IV <- simulate_mvn_p(mu_IV, sigma)

# V — SGP: single dominant gene
mu_V <- rep(0, m)
idx_V <- sample.int(m, 1)
mu_V[idx_V] <- 4.2
p_V <- simulate_mvn_p(mu_V, sigma)

# ---------------- build plotting data ----------------
df <- bind_rows(
  rank_curve(p_I,   "Archetype I — SDA"),
  rank_curve(p_II,  "Archetype II — CME"),
  rank_curve(p_III, "Archetype III — DPS"),
  rank_curve(p_IV,  "Archetype IV — HDS"),
  rank_curve(p_V,   "Archetype V — SGP")
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

k_1e3_sgp <- sum(sgp_dat$p < 1e-3)
k_gw_sgp  <- sum(sgp_dat$p < 2.5e-6)
top1 <- sgp_dat[1, , drop = FALSE]

p5 <- ggplot() +
  geom_line(data = sgp_dat,  aes(rank, mlog10p), linewidth = 0.7) +
  geom_point(data = sgp_dat, aes(rank, mlog10p), size = 0.55, alpha = 0.75) +
  geom_point(data = top1, aes(rank, mlog10p), size = 2) +
  coord_cartesian(ylim = c(0, 8.5)) +
  scale_x_continuous(labels = percent_format(accuracy = 1),
                     breaks = c(0, .25, .5, .75, 1)) +
  labs(
    title = "Archetype V — Single-Gene Proxy (SGP)",
    subtitle = "one dominant gene; remaining genes are mostly null",
    x = "Within-pathway gene rank (fraction)",
    y = expression(-log[10](p))
  ) +
  theme_natureish()
p5 <- add_thresholds(p5) +
  anno_box(0.12, 8.35,
           paste0("K(p<1e-3)=", k_1e3_sgp, "   K(p<2.5e-6)=", k_gw_sgp))

# ---------------- compose 2x3 layout for five archetypes ----------------
fig <- (p1 + p2 + p3) / (p4 + p5 + plot_spacer()) +
  plot_annotation(
    title = "Pathway signal archetypes (I-V) as gene-level p-value rank curves",
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
#ggsave("Fig_archetypes_layout.png", fig, width = 10.8, height = 7.0, units = "in", dpi = 600)

# ggsave("Fig_archetypes_rankcurves.pdf", fig, width = 7.2, height = 4.6, units = "in", device = cairo_pdf)
# ggsave("Fig_archetypes_rankcurves.png", fig, width = 7.2, height = 4.6, units = "in", dpi = 600)
