#!/usr/bin/env Rscript
# =============================================================================
# 13_spikein_eval.R  --  Evaluate staged spike-in recovery and calibration.
# =============================================================================
suppressMessages({
  library(data.table)
  library(ggplot2)
})

getenv1 <- function(k, default = "") {
  v <- Sys.getenv(k, default)
  if (!nzchar(v)) default else v
}

clean_gene_vec <- function(x) {
  x <- unique(as.character(x))
  x <- trimws(x)
  x[nzchar(x) & !is.na(x) & toupper(x) != "NA"]
}

calibration_table <- function(p, label) {
  p <- p[is.finite(p) & p > 0 & p <= 1]
  if (!length(p)) return(NULL)
  one <- function(a) {
    k <- sum(p <= a)
    ci <- binom.test(k, length(p), a)$conf.int
    data.table(
      null_set = label,
      alpha = a,
      n = length(p),
      observed = k / length(p),
      ci_lo = ci[1],
      ci_hi = ci[2],
      calibrated = (a >= ci[1] & a <= ci[2])
    )
  }
  rbind(one(0.05), one(0.01))
}

args <- commandArgs(trailingOnly = TRUE)
in_dir <- if (length(args) >= 1) args[1] else Sys.getenv("CATFISH_OUT")
pheno_dir <- if (length(args) >= 2) args[2] else Sys.getenv("PHENO_DIR")
out_dir <- if (length(args) >= 3) args[3] else Sys.getenv("SPIKE_SUMMARY_DIR", Sys.getenv("FINAL_DIR"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

truth <- fread(file.path(pheno_dir, "truth.tsv"))
spiked <- truth$spiked_pathway[1]
arch <- truth$archetype[1]
h2 <- truth$h2[1]
cat(sprintf("[eval] spiked=%s arch=%s h2=%.3f\n", spiked, arch, h2))
causal_genes <- clean_gene_vec(unlist(strsplit(getenv1("SPIKE_CAUSAL_GENES", truth$causal_genes[1]), ";", fixed = TRUE)))

files <- list.files(in_dir, pattern = "^pathway_pvals_perm_.*\\.csv$", full.names = TRUE)
if (!length(files)) stop(paste("no per-sim CSVs in", in_dir), call. = FALSE)
dat <- rbindlist(lapply(files, fread), fill = TRUE)
hp <- if ("omni_p_final" %in% names(dat)) "omni_p_final" else "omni_p_mvn"

per_sim <- dat[, {
  bh <- p.adjust(get(hp), method = "BH")
  bonf <- p.adjust(get(hp), method = "bonferroni")
  i <- which(pathway_id == spiked)
  data.table(
    spiked_p = if (length(i)) get(hp)[i] else NA_real_,
    spiked_fdr = if (length(i)) bh[i] else NA_real_,
    spiked_bonf = if (length(i)) bonf[i] else NA_real_,
    spiked_rank = if (length(i)) rank(get(hp))[i] else NA_real_,
    n_paths = .N
  )
}, by = perm]

power <- data.table(
  metric = c("nominal_p<0.05", "FDR<0.05", "Bonferroni<0.05"),
  power = c(
    mean(per_sim$spiked_p < 0.05, na.rm = TRUE),
    mean(per_sim$spiked_fdr < 0.05, na.rm = TRUE),
    mean(per_sim$spiked_bonf < 0.05, na.rm = TRUE)
  ),
  median_rank = median(per_sim$spiked_rank, na.rm = TRUE)
)
fwrite(power, file.path(out_dir, "power_summary.csv"))

# ---- overlap-aware null definitions -----------------------------------------
PATHWAY_FILE <- getenv1("PATHWAY_FILE")
pw_overlap <- NULL
if (nzchar(PATHWAY_FILE) && file.exists(PATHWAY_FILE)) {
  pw_raw <- fread(PATHWAY_FILE, check.names = FALSE)
  gene_cols <- intersect(c("Gene-name", "Gene-id"), names(pw_raw))
  if (length(gene_cols) && all(c("Pathway-id", "Pathway-name") %in% names(pw_raw))) {
    overlap_tbl <- rbindlist(lapply(gene_cols, function(gc) {
      pw <- unique(pw_raw[, .(
        pathway_id = `Pathway-id`,
        pathway_name = `Pathway-name`,
        gene_id = as.character(get(gc))
      )])
      pw <- pw[nzchar(gene_id) & !is.na(gene_id)]
      spiked_genes <- clean_gene_vec(pw[pathway_id == spiked, gene_id])
      if (!length(spiked_genes)) {
        data.table(
          gene_col = gc,
          n_spiked_genes = 0L,
          n_spiked_overlap = 0L,
          n_causal_genes = length(causal_genes),
          n_causal_overlap = 0L
        )
      } else {
        pathway_sets <- pw[, .(genes = list(unique(gene_id))), by = pathway_id]
        data.table(
          gene_col = gc,
          n_spiked_genes = length(spiked_genes),
          n_spiked_overlap = sum(vapply(pathway_sets$genes, function(g) any(g %chin% spiked_genes), logical(1))),
          n_causal_genes = length(causal_genes),
          n_causal_overlap = if (length(causal_genes)) {
            sum(vapply(pathway_sets$genes, function(g) any(g %chin% causal_genes), logical(1)))
          } else 0L
        )
      }
    }), fill = TRUE)

    overlap_tbl[, score := pmax(n_causal_overlap, n_spiked_overlap)]
    setorder(overlap_tbl, -score, -n_spiked_overlap, -n_causal_overlap)
    best_col <- overlap_tbl$gene_col[1]
    pw_overlap <- unique(pw_raw[, .(
      pathway_id = `Pathway-id`,
      pathway_name = `Pathway-name`,
      gene_id = as.character(get(best_col))
    )])
    pw_overlap <- pw_overlap[nzchar(gene_id) & !is.na(gene_id)]
    fwrite(overlap_tbl[, .(gene_col, n_spiked_genes, n_spiked_overlap, n_causal_genes, n_causal_overlap)],
           file.path(out_dir, "pathway_overlap_column_check.csv"))
    cat(sprintf("[eval] overlap gene column chosen: %s\n", best_col))
  }
}

null_sets <- list(
  naive_all_other = setdiff(unique(dat$pathway_id), spiked)
)

overlap_summary <- NULL
if (!is.null(pw_overlap)) {
  spiked_pathway_genes <- clean_gene_vec(pw_overlap[pathway_id == spiked, gene_id])
  pathway_sets <- pw_overlap[, .(genes = list(unique(gene_id))), by = pathway_id]

  if (length(spiked_pathway_genes)) {
    overlap_spiked <- pathway_sets[vapply(genes, function(g) any(g %chin% spiked_pathway_genes), logical(1)), pathway_id]
    null_sets[["nonoverlap_spiked_pathway"]] <- setdiff(unique(dat$pathway_id), overlap_spiked)
  }
  if (length(causal_genes)) {
    overlap_causal <- pathway_sets[vapply(genes, function(g) any(g %chin% causal_genes), logical(1)), pathway_id]
    null_sets[["nonoverlap_causal_genes"]] <- setdiff(unique(dat$pathway_id), overlap_causal)
  }

  overlap_summary <- rbindlist(lapply(names(null_sets), function(nm) {
    kept <- null_sets[[nm]]
    data.table(
      null_set = nm,
      n_pathways_kept = length(kept),
      n_rows_kept = sum(dat$pathway_id %chin% kept)
    )
  }), fill = TRUE)
  fwrite(overlap_summary, file.path(out_dir, "pathway_overlap_summary.csv"))
}

calib <- rbindlist(lapply(names(null_sets), function(nm) {
  kept <- null_sets[[nm]]
  calibration_table(dat[pathway_id %chin% kept][[hp]], nm)
}), fill = TRUE)
fwrite(calib, file.path(out_dir, "calibration_unspiked.csv"))

comp <- intersect(
  c("acat_p", "fisher_p", "minp_p", "stouffer_p_analytic", "tfisher_p", "tfisher_p_analytic"),
  names(dat)
)
if (length(comp)) {
  sp <- dat[pathway_id == spiked, ..comp]
  if (nrow(sp)) {
    winner <- comp[max.col(-as.matrix(sp), ties.method = "first")]
    wt <- as.data.table(table(component = winner))[order(-N)]
    wt[, pct := round(100 * N / sum(N), 1)]
    fwrite(wt, file.path(out_dir, "dominant_component.csv"))
  }
}

plot_null_set <- if ("nonoverlap_causal_genes" %in% calib$null_set) {
  "nonoverlap_causal_genes"
} else if ("nonoverlap_spiked_pathway" %in% calib$null_set) {
  "nonoverlap_spiked_pathway"
} else {
  "naive_all_other"
}
plot_null_paths <- null_sets[[plot_null_set]]
plot_null_p <- dat[pathway_id %chin% plot_null_paths][[hp]]
plot_null_p <- plot_null_p[is.finite(plot_null_p) & plot_null_p > 0 & plot_null_p <= 1]
if (!length(plot_null_p)) {
  plot_null_set <- "naive_all_other"
  plot_null_paths <- null_sets[[plot_null_set]]
  plot_null_p <- dat[pathway_id %chin% plot_null_paths][[hp]]
  plot_null_p <- plot_null_p[is.finite(plot_null_p) & plot_null_p > 0 & plot_null_p <= 1]
}
subtitle_bits <- calib[alpha == 0.05, .(null_set, observed)]
subtitle_txt <- paste(sprintf("%s=%.3f", subtitle_bits$null_set, subtitle_bits$observed), collapse = " | ")

plot_dt <- rbind(
  data.table(grp = "spiked pathway", val = -log10(per_sim$spiked_p)),
  data.table(grp = sprintf("null: %s", plot_null_set), val = -log10(sample(plot_null_p, min(length(plot_null_p), 5000))))
)
p <- ggplot(plot_dt, aes(val, fill = grp)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed", colour = "red") +
  labs(
    x = expression(-log[10](p)),
    fill = NULL,
    title = sprintf("Spike-in recovery: %s  (arch=%s, h2=%.2f)", spiked, arch, h2),
    subtitle = sprintf(
      "power(FDR<0.05)=%.2f | type-I@0.05: %s",
      power$power[power$metric == "FDR<0.05"],
      subtitle_txt
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "spikein_recovery.png"), p, width = 7, height = 4.5, dpi = 200)
ggsave(file.path(out_dir, "spikein_recovery.pdf"), p, width = 7, height = 4.5)
cat(sprintf("[eval] tables + figure -> %s\n", out_dir))
