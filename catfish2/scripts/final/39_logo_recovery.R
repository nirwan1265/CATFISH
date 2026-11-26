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
CTRL_LAB <- Sys.getenv("CONTROLS_LABEL", if (nzchar(Sys.getenv("CONTROLS_FILE"))) "Open Targets gold-standard genes" else "established genes")
if (nzchar(CF_)) {
  cf_ <- data.table::fread(CF_)
  if ("universe" %in% names(cf_)) cf_ <- cf_[cf_$universe == Sys.getenv("CONTROLS_UNIVERSE", "no_MHC"), ]
  tr_all <- unique(cf_$trait)
  CONTROLS <- lapply(setNames(tr_all, tr_all), function(t) unique(cf_$symbol[cf_$trait == t]))
  cat("controls from", CF_, ":", paste(names(CONTROLS), lengths(CONTROLS), collapse = ", "), "\n")
}
## SETS (comma-separated) picks the analyses: default the ten lipid sets; the
## benchmark traits are <trait>_eur (candidates in results/benchmark/candidates,
## WY results in results/benchmark/pathways). OUT_PREFIX names the outputs
## (default 39_logo_recovery; table tab_logo_recovery.tex, or tab_<OUT_PREFIX>.tex).
SETS <- strsplit(Sys.getenv("SETS", paste(c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_")), collapse = ",")), ",")[[1]]
OUT_PREFIX <- Sys.getenv("OUT_PREFIX", "39_logo_recovery")
TAB_NAME <- paste0("tab_", sub("^39_", "", OUT_PREFIX), ".tex")   # tab_logo_recovery.tex by default
BENCH <- c("cad", "brca", "cd", "uc", "af", "sle", "asthma")
## LOGO_CAP (default 1): a gene's pathway statistic without itself is never
## taken as stronger than with itself. Genes whose best pathway is the same share
## one tied pathway rank (averaged); without the cap, a gene whose removal leaves
## its pathway marginally stronger (it diluted the pathway) would jump from the
## middle of that tie to its head, a gain from tie-breaking, not from evidence.
## LOGO_CAP=0 reproduces the uncapped version (until 2026-11-25).
LOGO_CAP <- Sys.getenv("LOGO_CAP", "1") == "1"
## UNIVERSE_FILE (optional): rank only within the genes of each trait listed in
## this file (columns trait, entrez), e.g. results/45_universe_no_MHC.csv from
## final/45_benchmark_topk_all.R: the genes scored by both CATFISH and PoPS,
## outside the MHC, as in Table tab:benchmark_topk_all. Scores are not
## recomputed (as 45): each layer's ranking is re-ranked within those genes.
## Pathway evidence of a gene left out is still ranked among all genes'
## pathway statistics, as in 06. Unset: all genes of the 06 file (until 2026-11-26).
UNI_F <- Sys.getenv("UNIVERSE_FILE", "")
UNI <- if (nzchar(UNI_F)) { u_ <- fread(UNI_F, colClasses = list(character = "entrez")); split(u_$entrez, u_$trait) } else NULL
if (!is.null(UNI)) cat("universe from", UNI_F, ":", paste(names(UNI), lengths(UNI), collapse = ", "), "\n")

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
  f06 <- if (tr %in% BENCH) file.path(RES, "benchmark", "candidates", sprintf("06_candidate_genes_%s.csv", tr)) else
    file.path(RES, sprintf("06_candidate_genes_%s.csv", tag))
  g <- fread(f06, colClasses = list(character = "gene"))
  g[, `:=`(S_GWAS = z0(S_GWAS), S_MAGMA = z0(S_MAGMA), S_PATH = z0(S_PATH), B_hub = z0(B_hub))]
  g[, r_gm   := avg_rank_desc(S_GWAS + S_MAGMA)]
  g[, r_gmh  := avg_rank_desc(S_GWAS + S_MAGMA + W_HUB * B_hub)]
  g[, r_full := avg_rank_desc(S_GWAS + S_MAGMA + W_PATH * S_PATH + W_HUB * B_hub)]
  score <- g$S_GWAS + g$S_MAGMA + W_PATH * g$S_PATH + W_HUB * g$B_hub
  if (max(abs(score - g$Score)) > 1e-8) stop(A, ": recomputed score does not reproduce 06's Score", call. = FALSE)
  pstat <- suppressWarnings(as.numeric(g$p_pathway)); N_PATH <- sum(!is.na(pstat))
  ## ranking universe: all genes, or those of UNIVERSE_FILE (layer ranks re-ranked within it)
  if (is.null(UNI)) {
    g[, inU := TRUE]; g[, `:=`(k_gwas = rank_gwas, k_magma = rank_magma, k_full = rank)]
  } else {
    if (is.null(UNI[[tr]])) stop(A, ": trait ", tr, " not in ", UNI_F, call. = FALSE)
    g[, inU := gene %in% UNI[[tr]]]
    g[, `:=`(k_gwas = NA_real_, k_magma = NA_real_, k_full = NA_real_)]
    g[inU == TRUE, `:=`(k_gwas = frank(rank_gwas, ties.method = "average", na.last = TRUE),
                        k_magma = frank(rank_magma, ties.method = "average", na.last = TRUE),
                        k_full = avg_rank_desc(Score))]
    cat(sprintf("%s: %d of %d genes in the universe\n", A, sum(g$inU), nrow(g)))
  }
  g[, `:=`(r_gm = NA_real_, r_gmh = NA_real_)]
  g[inU == TRUE, r_gm  := avg_rank_desc(S_GWAS + S_MAGMA)]
  g[inU == TRUE, r_gmh := avg_rank_desc(S_GWAS + S_MAGMA + W_HUB * B_hub)]

  wy <- fread(file.path(if (tr %in% BENCH) file.path(RES, "benchmark", "pathways") else WY,
                        sprintf("%s_reactome_5comp_westfall_young.csv", A)))
  zp <- setNames(wy$z_std, wy$pathway_id)
  fl <- file.path(LOGO, sprintf("%s_logo.csv", A))
  lo <- if (file.exists(fl)) fread(fl, colClasses = list(character = "removed_gene")) else NULL
  if (!is.null(lo)) lo[!nzchar(removed_gene), removed_gene := NA_character_]
  if (is.null(lo)) cat(sprintf("%s: no LOGO file (%s); LOGO column left empty\n", A, fl))

  for (s in CONTROLS[[tr]]) {
    i <- match(s, g$symbol); if (is.na(i) || !g$inU[i]) next
    id <- g$gene[i]
    r_logo <- NA_real_; best_logo <- NA_character_; mc <- NA_real_; p_logo <- NA_real_
    if (!is.null(lo)) {
      if (is.na(pstat[i])) { r_logo <- g$k_full[i]          # in no pathway: nothing to remove
      } else {
        m <- lo[removed_gene == id & !is.na(z_std)]
        base <- lo[is.na(removed_gene), .(base_pathway = pathway_id, z_base = z_std)]
        m <- merge(m, base, by = "base_pathway")
        if (nrow(m)) {
          m[, z_logo := zp[base_pathway] + (z_std - z_base)]
          mc <- max(abs(zp[m$base_pathway] - m$z_base))  # Monte Carlo difference from the paper's run
          j <- which.min(m$z_logo); p_logo <- m$z_logo[j]; best_logo <- m$base_pathway[j]
          if (LOGO_CAP) p_logo <- max(p_logo, pstat[i])   # no gain from its own removal
          rk <- rank_in(p_logo, pstat[-i][!is.na(pstat[-i])], low = TRUE)
          sp <- -log10(rk / N_PATH)
          sc <- g$S_GWAS[i] + g$S_MAGMA[i] + W_PATH * sp + W_HUB * g$B_hub[i]
          r_logo <- rank_in(sc, score[g$inU & seq_along(score) != i])
        } else r_logo <- NA_real_   # its pathways all fell below the size limit without it
      }
    }
    genes_out[[length(genes_out) + 1L]] <- data.table(
      set = A, gene = s, rank_gwas = g$k_gwas[i], rank_magma = g$k_magma[i],
      rank_gwas_magma = g$r_gm[i], rank_gwas_magma_hub = g$r_gmh[i], rank_catfish = g$k_full[i],
      rank_catfish_logo = r_logo, best_pathway = g$best_pathway[i], z_best = pstat[i],
      best_pathway_logo = best_logo, z_best_logo = p_logo, mc_diff = mc)
  }
}
G <- rbindlist(genes_out)
fwrite(G, file.path(RES, paste0(OUT_PREFIX, "_genes.csv")))

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
fwrite(S, file.path(RES, paste0(OUT_PREFIX, "_summary.csv")))

## ---- table ------------------------------------------------------------------------
TL <- c(ldl = "LDL-C", hdl = "HDL-C", cad = "CAD", brca = "Breast cancer", cd = "Crohn's disease",
        uc = "Ulcerative colitis", af = "Atrial fibrillation", sle = "SLE", asthma = "Asthma")
lab <- function(A) sprintf("%s & %s", TL[[sub("_.*", "", A)]],
                           toupper(sub(".*_", "", A)))
fm <- function(m, t) if (is.na(m)) "--" else sprintf("%s (%d)", trimws(formatC(m, format = "fg", digits = 4)), t)
WHERE <- if (all(sub("_.*", "", SETS) %in% BENCH)) sprintf(", in the %s further European GWAS", c("one", "two", "three", "four", "five", "six", "seven", "eight", "nine")[length(SETS)]) else ", in the lipid analyses"
L <- c("\\begin{table*}[!t]", "\\centering",
  paste0("\\caption{Where the gain of the gene score comes from", WHERE, ". Median rank of the ", CTRL_LAB, " ",
         "among the ranked genes (number in the top 20 in parentheses) when genes are ranked by the ",
         "GWAS layer alone, the MAGMA layer alone, the two together (pathway weight $0$), the full ",
         "CATFISH score, and the full score with each gene's pathway evidence recomputed ",
         "without that gene (leave one gene out; never stronger than with it).",
         if (!is.null(UNI)) " Genes are ranked within those scored by both CATFISH and PoPS, outside the MHC, as in Table~\\ref{tab:benchmark_topk_all}." else "", "}"),
  sprintf("\\label{tab:%s}", sub("\\.tex$", "", sub("^tab_", "", TAB_NAME))), "\\small",
  "\\begin{tabular}{llrrrrr}", "\\hline",
  paste0("\\textbf{Trait} & \\textbf{Group} & GWAS & MAGMA & GWAS + MAGMA & CATFISH & ",
         "\\makecell[r]{CATFISH,\\\\leave one gene out} \\\\"), "\\hline")
for (k in seq_len(nrow(S))) { r <- S[k]
  L <- c(L, sprintf("%s & %s & %s & %s & %s & %s \\\\", lab(r$set), fm(r$gwas, r$top20_gwas),
                    fm(r$magma, r$top20_magma), fm(r$gwas_magma, r$top20_gwas_magma),
                    fm(r$catfish, r$top20_catfish), fm(r$catfish_logo, r$top20_catfish_logo))) }
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, TAB_NAME))

options(width = 180)
print(S[, .(set, gwas, magma, gwas_magma, gwas_magma_hub, catfish, catfish_logo,
            top20_gwas_magma, top20_catfish, top20_catfish_logo, max_mc_diff)], row.names = FALSE)
cat("\nwrote", file.path(RES, paste0(OUT_PREFIX, "_summary.csv")), "and", file.path(AUT, TAB_NAME), "\n")
