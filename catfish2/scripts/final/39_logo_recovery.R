#!/usr/bin/env Rscript
## =============================================================================
## 39_logo_recovery.R -- is the gain of the gene score over GWAS and MAGMA
## pathway-informed? Two checks on the positive controls of the ten lipid sets:
##
##  1. Layers. The score of Section 2.5 recomputed from 06's own columns with
##     the pathway weight set to 0 (GWAS + MAGMA), with the hub term only, and
##     in full (pathway 0.5, hub 0.3; reproduces 06's rank exactly).
##  2. Leave one gene out (LOGO). Each control is given the pathway evidence of
##     its pathways with itself removed: for every tested pathway P containing
##     it, z_std(P) + [z(P without g) - z(P)], the bracket from one WY run on
##     the same null draws (scripts/hpc/gc_real/wy_logo.R); its best pathway is
##     the strongest of these. Its pathway rank among all genes, score and
##     overall rank are recomputed with every other gene unchanged. If the
##     control still ranks high, the lift comes from the other genes of its
##     pathways, not from its own association.
##
## Reads  catfish2/results/06_candidate_genes_<tag>.csv
##        <BASE>/results/human_LDL/catfish/<A>_reactome_5comp_westfall_young.csv
##        catfish2/results/gc_real_logo/<A>_logo.csv   (rsync of gc_real/out_logo)
## Writes catfish2/results/39_logo_recovery_genes.csv, 39_logo_recovery_summary.csv
##        catfish2/tex/auto/tab_logo_recovery.tex
## Usage (from catfish2/scripts/final):  Rscript 39_logo_recovery.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- Sys.getenv("RES_DIR",  file.path(BASE, "catfish2/results"))
WY   <- Sys.getenv("WY_DIR",   file.path(BASE, "results/human_LDL/catfish"))
LOGO <- Sys.getenv("LOGO_DIR", file.path(RES, "gc_real_logo"))
AUT  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
W_PATH <- as.numeric(Sys.getenv("W_PATH", "1")); W_HUB <- as.numeric(Sys.getenv("W_HUB", "0"))   # as 06_candidate_genes.R (until 2026-11-25: 0.5, 0.3)

CONTROLS <- list(
  ldl = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
          "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
  hdl = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
          "PLTP","LPL","ANGPTL4","CUBN"))

## CONTROLS_FILE (optional): evaluate on another gene list instead of the
## pre-specified one, e.g. the Open Targets gold-standard genes written by
## final/45_benchmark_topk_all.R (results/45_gold_genes.csv: trait, universe,
## ENSGID, symbol, entrez; rows of CONTROLS_UNIVERSE, default no_MHC, are used).
## CONTROLS_LABEL names the list in figure and table text.
CF_ <- Sys.getenv("CONTROLS_FILE", "")
CTRL_LAB <- Sys.getenv("CONTROLS_LABEL", "established genes")
if (nzchar(CF_)) {
  cf_ <- data.table::fread(CF_)
  if ("universe" %in% names(cf_)) cf_ <- cf_[cf_$universe == Sys.getenv("CONTROLS_UNIVERSE", "no_MHC"), ]
  CONTROLS <- lapply(setNames(names(CONTROLS), names(CONTROLS)), function(t) unique(cf_$symbol[cf_$trait == t]))
  cat("controls from", CF_, ":", paste(names(CONTROLS), lengths(CONTROLS), collapse = ", "), "\n")
}
SETS <- c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_"))

z0 <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }
## rank of value v among `others` (larger = better unless low = TRUE); average ties
rank_in <- function(v, others, low = FALSE) {
  if (low) { v <- -v; others <- -others }
  1 + sum(others > v) + 0.5 * sum(others == v)
}
avg_rank_desc <- function(x) frank(-x, ties.method = "average")

genes_out <- list(); summ <- list()
for (A in SETS) {
  tr <- sub("_.*", "", A); anc <- sub(".*_", "", A); tag <- if (anc == "eur") tr else A
  g <- fread(file.path(RES, sprintf("06_candidate_genes_%s.csv", tag)), colClasses = list(character = "gene"))
  g[, `:=`(S_GWAS = z0(S_GWAS), S_MAGMA = z0(S_MAGMA), S_PATH = z0(S_PATH), B_hub = z0(B_hub))]
  g[, r_gm   := avg_rank_desc(S_GWAS + S_MAGMA)]
  g[, r_gmh  := avg_rank_desc(S_GWAS + S_MAGMA + W_HUB * B_hub)]
  g[, r_full := avg_rank_desc(S_GWAS + S_MAGMA + W_PATH * S_PATH + W_HUB * B_hub)]
  score <- g$S_GWAS + g$S_MAGMA + W_PATH * g$S_PATH + W_HUB * g$B_hub
  if (max(abs(score - g$Score)) > 1e-8) stop(A, ": recomputed score does not reproduce 06's Score", call. = FALSE)
  pstat <- suppressWarnings(as.numeric(g$p_pathway)); N_PATH <- sum(!is.na(pstat))

  wy <- fread(file.path(WY, sprintf("%s_reactome_5comp_westfall_young.csv", A)))
  zp <- setNames(wy$z_std, wy$pathway_id)
  fl <- file.path(LOGO, sprintf("%s_logo.csv", A))
  lo <- if (file.exists(fl)) fread(fl, colClasses = list(character = "removed_gene")) else NULL
  if (!is.null(lo)) lo[!nzchar(removed_gene), removed_gene := NA_character_]
  if (is.null(lo)) cat(sprintf("%s: no LOGO file (%s); LOGO column left empty\n", A, fl))

  for (s in CONTROLS[[tr]]) {
    i <- match(s, g$symbol); if (is.na(i)) next
    id <- g$gene[i]
    r_logo <- NA_real_; best_logo <- NA_character_; mc <- NA_real_; p_logo <- NA_real_
    if (!is.null(lo)) {
      if (is.na(pstat[i])) { r_logo <- g$rank[i]          # in no pathway: nothing to remove
      } else {
        m <- lo[removed_gene == id & !is.na(z_std)]
        base <- lo[is.na(removed_gene), .(base_pathway = pathway_id, z_base = z_std)]
        m <- merge(m, base, by = "base_pathway")
        if (nrow(m)) {
          m[, z_logo := zp[base_pathway] + (z_std - z_base)]
          mc <- max(abs(zp[m$base_pathway] - m$z_base))  # Monte Carlo difference from the paper's run
          j <- which.min(m$z_logo); p_logo <- m$z_logo[j]; best_logo <- m$base_pathway[j]
          rk <- rank_in(p_logo, pstat[-i][!is.na(pstat[-i])], low = TRUE)
          sp <- -log10(rk / N_PATH)
          sc <- g$S_GWAS[i] + g$S_MAGMA[i] + W_PATH * sp + W_HUB * g$B_hub[i]
          r_logo <- rank_in(sc, score[-i])
        } else r_logo <- NA_real_   # its pathways all fell below the size limit without it
      }
    }
    genes_out[[length(genes_out) + 1L]] <- data.table(
      set = A, gene = s, rank_gwas = g$rank_gwas[i], rank_magma = g$rank_magma[i],
      rank_gwas_magma = g$r_gm[i], rank_gwas_magma_hub = g$r_gmh[i], rank_catfish = g$rank[i],
      rank_catfish_logo = r_logo, best_pathway = g$best_pathway[i], z_best = pstat[i],
      best_pathway_logo = best_logo, z_best_logo = p_logo, mc_diff = mc)
  }
}
G <- rbindlist(genes_out)
fwrite(G, file.path(RES, "39_logo_recovery_genes.csv"))

med <- function(x) as.numeric(stats::median(x, na.rm = TRUE))
t20 <- function(x) sum(x <= 20, na.rm = TRUE)
S <- G[, .(n = .N,
           gwas = med(rank_gwas), magma = med(rank_magma), gwas_magma = med(rank_gwas_magma),
           gwas_magma_hub = med(rank_gwas_magma_hub), catfish = med(rank_catfish),
           catfish_logo = med(rank_catfish_logo),
           top20_gwas = t20(rank_gwas), top20_magma = t20(rank_magma),
           top20_gwas_magma = t20(rank_gwas_magma), top20_catfish = t20(rank_catfish),
           top20_catfish_logo = t20(rank_catfish_logo),
           max_mc_diff = if (all(is.na(mc_diff))) NA_real_ else max(mc_diff, na.rm = TRUE)), by = set]
fwrite(S, file.path(RES, "39_logo_recovery_summary.csv"))

## ---- table ------------------------------------------------------------------------
lab <- function(A) sprintf("%s & %s", c(ldl = "LDL-C", hdl = "HDL-C")[[sub("_.*", "", A)]],
                           toupper(sub(".*_", "", A)))
fm <- function(m, t) if (is.na(m)) "--" else sprintf("%s (%d)", trimws(formatC(m, format = "fg", digits = 4)), t)
L <- c("\\begin{table*}[!t]", "\\centering",
  paste0("\\caption{Where the gain of the gene score comes from. Median rank of the established ",
         "genes among all genes (number in the top 20 in parentheses) when genes are ranked by the ",
         "GWAS layer alone, the MAGMA layer alone, the two together (pathway weight $0$), the full ",
         "CATFISH score, and the full score with each established gene's pathway evidence recomputed ",
         "without that gene (leave one gene out).}"),
  "\\label{tab:logo_recovery}", "\\small",
  "\\begin{tabular}{llrrrrr}", "\\hline",
  paste0("\\textbf{Trait} & \\textbf{Group} & GWAS & MAGMA & GWAS + MAGMA & CATFISH & ",
         "\\makecell[r]{CATFISH,\\\\leave one gene out} \\\\"), "\\hline")
for (k in seq_len(nrow(S))) { r <- S[k]
  L <- c(L, sprintf("%s & %s & %s & %s & %s & %s \\\\", lab(r$set), fm(r$gwas, r$top20_gwas),
                    fm(r$magma, r$top20_magma), fm(r$gwas_magma, r$top20_gwas_magma),
                    fm(r$catfish, r$top20_catfish), fm(r$catfish_logo, r$top20_catfish_logo))) }
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, "tab_logo_recovery.tex"))

options(width = 180)
print(S[, .(set, gwas, magma, gwas_magma, gwas_magma_hub, catfish, catfish_logo,
            top20_gwas_magma, top20_catfish, top20_catfish_logo, max_mc_diff)], row.names = FALSE)
cat("\nwrote", file.path(RES, "39_logo_recovery_summary.csv"), "and", file.path(AUT, "tab_logo_recovery.tex"), "\n")
