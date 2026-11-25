#!/usr/bin/env Rscript
## =============================================================================
## 47_bap_logo.R -- leave one gene out (LOGO) for the BAP stem-composition gene
## score (Section 3.10), as 39_logo_recovery.R for the lipid controls. Each gene
## checked by hpc/bap_gwas/06_logo_bap.R (the pathway-supported genes of each
## trait's top 100, SbALDH3E1 and the thioredoxin reductase among them) is given
## the pathway evidence of its pathways with itself removed: for every tested
## pathway P containing it, z_std(P) + [z(P without g) - z(P)], the bracket from
## one WY run on the same null draws; its pathway statistic is the strongest of
## these. Its pathway score, gene score and rank (among all genes, and among the
## genes with evidence in all three layers) are recomputed with every other gene
## unchanged. If a gene still ranks high, its lift comes from the other genes of
## its pathways, not from its own association.
## Reads  catfish2/results/bap_traits/prioritize/<trait>_gene_prioritize.csv
##        catfish2/results/bap_traits/<trait>/<trait>_wy.csv
##        catfish2/results/bap_traits/logo/<trait>_logo.csv  (rsync of BAP_gwas/logo)
## Writes catfish2/results/bap_traits/logo/47_bap_logo_genes.csv, 47_bap_logo_summary.csv
##        catfish2/tex/auto/bap_logo.tex (+ bap_logo_stub.tex)
##   Rscript 47_bap_logo.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- file.path(BASE, "catfish2/results/bap_traits")
TEX  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
TR   <- c(Brix = "Brix_maturity", NDF = "NDF_pct_DM", NFC = "NFC_pct_DM")
W_PATH <- as.numeric(Sys.getenv("W_PATH", "1")); W_HUB <- as.numeric(Sys.getenv("W_HUB", "0"))   # gene_prioritize() defaults (until 2026-11-25: 0.5, 0.3)
ALDH <- "SORBI_3004G300800"; TRXR <- "SORBI_3004G300900"
LOGO_CAP <- Sys.getenv("LOGO_CAP", "1") == "1"   # as 39_logo_recovery.R: without itself never stronger than with itself (0 = uncapped, until 2026-11-25)
rank_in <- function(v, others, low = FALSE) {      # rank of v among others, average ties
  if (low) { v <- -v; others <- -others }
  1 + sum(others > v) + 0.5 * sum(others == v)
}
z0 <- function(x) { x <- suppressWarnings(as.numeric(x)); x[is.na(x)] <- 0; x }

G <- list(); S <- list(); M <- list()
for (k in names(TR)) {
  t <- TR[[k]]
  p  <- fread(file.path(RES, "prioritize", sprintf("%s_gene_prioritize.csv", t)))
  wy <- fread(file.path(RES, t, sprintf("%s_wy.csv", t)))
  fl <- file.path(RES, "logo", sprintf("%s_logo.csv", t))
  if (!file.exists(fl)) { cat("missing", fl, "\n"); next }
  lo <- fread(fl, colClasses = list(character = "removed_gene"))
  lo[!nzchar(removed_gene), removed_gene := NA_character_]
  p[, `:=`(gwas_score = z0(gwas_score), magma_score = z0(magma_score),
           pathway_score = z0(pathway_score), hub_bonus = z0(hub_bonus))]
  score <- p$gwas_score + p$magma_score + W_PATH * p$pathway_score + W_HUB * p$hub_bonus
  if (max(abs(score - p$score)) > 1e-8) stop(t, ": recomputed score does not reproduce gene_prioritize()", call. = FALSE)
  pstat <- suppressWarnings(as.numeric(p$pathway_stat)); N_PATH <- sum(!is.na(pstat))
  three <- p$n_layers == 3
  zp <- setNames(wy$z_std, wy$pathway_id)
  base <- lo[is.na(removed_gene), .(base_pathway = pathway_id, z_base = z_std)]
  gm <- frank(-(p$gwas_score + p$magma_score), ties.method = "average")
  for (id in unique(na.omit(lo$removed_gene))) {
    i <- match(id, p$gene); if (is.na(i)) next
    m <- merge(lo[removed_gene == id & !is.na(z_std)], base, by = "base_pathway")
    r_logo <- NA_real_; r3_logo <- NA_real_; best <- NA_character_
    if (nrow(m)) {
      m[, z_logo := zp[base_pathway] + (z_std - z_base)]
      j <- which.min(m$z_logo); best <- m$base_pathway[j]
      zl <- if (LOGO_CAP && !is.na(pstat[i])) max(m$z_logo[j], pstat[i]) else m$z_logo[j]
      sp <- -log10(rank_in(zl, pstat[-i][!is.na(pstat[-i])], low = TRUE) / N_PATH)
      sc <- p$gwas_score[i] + p$magma_score[i] + W_PATH * sp + W_HUB * p$hub_bonus[i]
      r_logo  <- rank_in(sc, score[-i])
      r3_logo <- rank_in(sc, score[-i][three[-i]])
    }
    G[[length(G) + 1L]] <- data.table(trait = t, gene = id, chr = p$chr[i], start = p$start[i],
      best_pathway = p$best_pathway_name[i], rank_gwas = p$gwas_rank[i], rank_magma = p$magma_rank[i],
      rank_gwas_magma = gm[i], rank_catfish = p$rank[i], rank_catfish_logo = r_logo,
      rank3_catfish = rank_in(score[i], score[-i][three[-i]]), rank3_catfish_logo = r3_logo,
      best_pathway_logo = best)
  }
}
X <- rbindlist(G)
dir.create(file.path(RES, "logo"), showWarnings = FALSE)
fwrite(X, file.path(RES, "logo", "47_bap_logo_genes.csv"))
S <- X[, .(n_genes = .N, median_rank_gwas_magma = as.numeric(median(rank_gwas_magma)), median_rank_catfish = as.numeric(median(rank_catfish)),
           median_rank_logo = as.numeric(median(rank_catfish_logo, na.rm = TRUE)),
           top100_catfish = sum(rank_catfish <= 100), top100_logo = sum(rank_catfish_logo <= 100, na.rm = TRUE)),
       by = trait]
fwrite(S, file.path(RES, "logo", "47_bap_logo_summary.csv"))
options(width = 200); print(S, row.names = FALSE)
cat("\nSbALDH3E1 and the thioredoxin reductase:\n")
print(X[gene %in% c(ALDH, TRXR), .(trait, gene, rank_gwas_magma, rank_catfish, rank_catfish_logo,
                                   rank3_catfish, rank3_catfish_logo)], row.names = FALSE)

## ---- macros ------------------------------------------------------------------------------
f0 <- function(v) if (is.na(v)) "--" else if (abs(v - round(v)) < 1e-9) format(round(v), big.mark = "{,}") else sprintf("%.1f", v)
for (k in names(TR)) {
  for (g in c(ALDH, TRXR)) {
    tg <- if (g == ALDH) "Aldh" else "Trx"; r <- X[trait == TR[[k]] & gene == g]
    if (!nrow(r)) next
    M[[sprintf("nBapLogo%s%s", tg, k)]]   <- f0(r$rank_catfish_logo)
    M[[sprintf("nBapLogoTl%s%s", tg, k)]] <- f0(r$rank3_catfish_logo)
    M[[sprintf("nBapGm%s%s", tg, k)]]     <- f0(r$rank_gwas_magma)
  }
  s <- S[trait == TR[[k]]]
  if (nrow(s)) {
    M[[sprintf("nBapLogoN%s", k)]] <- s$n_genes
    M[[sprintf("nBapLogoKeep%s", k)]] <- s$top100_logo
    M[[sprintf("nBapLogoMedGm%s", k)]] <- f0(s$median_rank_gwas_magma)
    M[[sprintf("nBapLogoMedCat%s", k)]] <- f0(s$median_rank_catfish)
    M[[sprintf("nBapLogoMed%s", k)]] <- f0(s$median_rank_logo)
  }
}
ALL <- c(as.vector(outer(c("nBapLogoAldh", "nBapLogoTlAldh", "nBapGmAldh", "nBapLogoTrx", "nBapLogoTlTrx", "nBapGmTrx"),
                         names(TR), paste0)),
         as.vector(outer(c("nBapLogoN", "nBapLogoKeep", "nBapLogoMedGm", "nBapLogoMedCat", "nBapLogoMed"), names(TR), paste0)))
writeLines(c("%% Placeholders for the macros written by 47_bap_logo.R",
             sprintf("\\providecommand{\\%s}{\\textbf{??}}", ALL)), file.path(TEX, "bap_logo_stub.tex"))
writeLines(c("%% written by 47_bap_logo.R -- do not edit",
             sprintf("\\renewcommand{\\%s}{%s}", names(M), unlist(M))), file.path(TEX, "bap_logo.tex"))
cat("wrote", file.path(RES, "logo", "47_bap_logo_genes.csv"), "and", file.path(TEX, "bap_logo.tex"), "\n")
