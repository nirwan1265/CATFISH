#!/usr/bin/env Rscript
## =============================================================================
## 44_bap_chr4_numbers.R -- numbers and table for the sorghum stem-composition
## section (BAP Brix at maturity, NDF, NFC; GEMMA with 3 PCs -> MAGMA -> CATFISH
## -> gene score), and the chr4 63.9-64.0 Mb locus.
## Reads   catfish2/data/bap/BAP_pheno_vcf2gwas.csv
##         catfish2/results/bap_traits/<trait>/<trait>_wy.csv        (hpc/bap_gwas/03)
##         catfish2/results/bap_traits/prioritize/<trait>_gene_prioritize.csv (hpc/bap_gwas/04)
##         catfish2/results/bap_traits/chr4_locus/                   (hpc/bap_gwas/05, final/43)
## Writes  catfish2/tex/auto/bap_chr4_stub.tex  (\newcommand placeholders)
##         catfish2/tex/auto/bap_chr4.tex       (\renewcommand values; a missing
##                                               input leaves its macros at ??)
##         catfish2/tex/auto/tab_bap_chr4.tex   (genes of the locus)
##   Rscript 44_bap_chr4_numbers.R
## Env: CATFISH_BASE, MAGMA_WINDOW (25000), LOCUS_FROM / LOCUS_TO (63900000 / 64060000)
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
RES  <- file.path(BASE, "catfish2/results/bap_traits")
LOC  <- file.path(RES, "chr4_locus"); PRD <- file.path(RES, "prioritize")
PHE  <- file.path(BASE, "catfish2/data/bap/BAP_pheno_vcf2gwas.csv")
TEX  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
MW   <- as.numeric(Sys.getenv("MAGMA_WINDOW", "25000"))
LFROM <- as.numeric(Sys.getenv("LOCUS_FROM", "63900000")); LTO <- as.numeric(Sys.getenv("LOCUS_TO", "64060000"))
TR   <- c(Brix = "Brix_maturity", NDF = "NDF_pct_DM", NFC = "NFC_pct_DM")
ALDH <- "SORBI_3004G300800"; TRXR <- "SORBI_3004G300900"
dir.create(TEX, showWarnings = FALSE, recursive = TRUE)

M <- list()                                    # macro name -> value (character)
put <- function(k, v) M[[k]] <<- as.character(v)
big <- function(x) gsub(",", "{,}", formatC(x, format = "d", big.mark = ","))
rk  <- function(x) big(floor(as.numeric(x) + 0.5))      # tied (mid) ranks rounded up
pf  <- function(p, d = 2) {                    # p-value in LaTeX math, without $
  if (!is.finite(p)) return("??")
  if (p >= 0.001) return(formatC(p, digits = d, format = "fg", flag = "#") |> sub("\\.$", "", x = _))
  e <- floor(log10(p)); m <- p / 10^e
  if (round(m, 1) >= 10) { m <- m / 10; e <- e + 1 }
  sprintf("%.1f\\times10^{%d}", m, e)
}
f2  <- function(x, d = 2) formatC(x, digits = d, format = "f")
tex_esc <- function(s) gsub("&", "\\\\&", gsub("%", "\\\\%", s))
pw_name <- function(s) { s <- gsub("&beta;", "$\\\\beta$", s); s <- gsub("&alpha;", "$\\\\alpha$", s)
  s <- gsub("<i>|</i>", "", s); gsub("<[^>]+>", "", s) }

## ---- phenotypes -------------------------------------------------------------------
ph <- fread(PHE); setnames(ph, 1, "id")
raw <- if (file.exists(file.path(LOC, "leads_geno.raw"))) fread(file.path(LOC, "leads_geno.raw")) else NULL
if (!is.null(raw)) {
  gid <- raw$IID
  put("nBapNlines", nrow(raw))
  for (k in names(TR)) put(sprintf("nBapN%s", k), sum(ph$id %in% gid & is.finite(ph[[TR[k]]])))
}
cr <- function(a, b) cor(ph[[a]], ph[[b]], use = "pairwise.complete.obs")
put("nBapRBrixNDF", f2(cr(TR["Brix"], TR["NDF"]))); put("nBapRBrixNFC", f2(cr(TR["Brix"], TR["NFC"])))
put("nBapRNDFNFC", f2(cr(TR["NDF"], TR["NFC"])))
put("nBapRBrixADF", f2(cr(TR["Brix"], "ADF_pct_DM"))); put("nBapRBrixLig", f2(cr(TR["Brix"], "Lignin_pct_DM")))

## ---- pathways ---------------------------------------------------------------------
sig_list <- function(w) {
  s <- w[p_fwer_std <= 0.05][order(p_fwer_std)]
  if (!nrow(s)) return("none")
  paste(sprintf("%s ($p=%s$)", tex_esc(pw_name(s$pathway_name)), vapply(s$p_fwer_std, pf, "", d = 2)), collapse = ", ")
}
for (k in names(TR)) {
  f <- file.path(RES, TR[k], sprintf("%s_wy.csv", TR[k]))
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  w <- fread(f)
  put(sprintf("nBapGC%s", k), f2(w$gc_shift[1]))
  put(sprintf("nBapSig%s", k), sum(w$p_fwer_std <= 0.05))
  put(sprintf("nBapComp%s", k), sum(w$p_fwer_comp <= 0.05, na.rm = TRUE))
  put(sprintf("nBapSigList%s", k), sig_list(w))
  put("nBapNpw", nrow(w))
  ba <- w[pathway_id == "PWY-3981"]; tx <- w[pathway_id == "THIOREDOX-PWY"]
  if (nrow(ba)) { put(sprintf("nBapBala%s", k), pf(ba$p_fwer_std)); put(sprintf("nBapBalaComp%s", k), pf(ba$p_fwer_comp))
                  put("nBapBalaN", ba$n_genes) }
  if (nrow(tx)) { put(sprintf("nBapTrxPw%s", k), pf(tx$p_fwer_std)); put(sprintf("nBapTrxPwComp%s", k), pf(tx$p_fwer_comp)) }
}

## ---- gene score -------------------------------------------------------------------
PR <- list()
for (k in names(TR)) {
  f <- file.path(PRD, sprintf("%s_gene_prioritize.csv", TR[k]))
  if (!file.exists(f)) { cat("missing", f, "\n"); next }
  x <- fread(f); x[, pos := start + MW]; PR[[k]] <- x
  three <- x[n_layers == 3][order(-score)][, tl := .I]
  put("nBapNgenes", big(nrow(x)))
  put(sprintf("nBapNthree%s", k), big(nrow(three)))
  for (g in c(ALDH, TRXR)) {
    tag <- if (g == ALDH) "Aldh" else "Trx"
    put(sprintf("nBap%s%s", tag, k), rk(x[gene == g, rank]))
    put(sprintf("nBap%sTl%s", tag, k), rk(three[gene == g, tl]))
    put(sprintf("nBap%sGwas%s", tag, k), rk(x[gene == g, gwas_rank]))
    put(sprintf("nBap%sMagmaP%s", tag, k), pf(x[gene == g, magma_p]))
  }
  top <- x[which.max(gwas_neglog10p)]
  put(sprintf("nBapTopLoc%s", k), sprintf("chromosome~%s at $%.1f$~Mb", top$chr, top$pos / 1e6))
  put(sprintf("nBapTopP%s", k), pf(10^-top$gwas_neglog10p))
  blk <- x[as.character(chr) == "4" & pos >= LFROM & pos <= LTO]
  put(sprintf("nBapLocP%s", k), pf(10^-max(blk$gwas_neglog10p)))
}
if (length(PR) == 3) {
  W <- Reduce(function(a, b) merge(a, b, by = "gene"),
              lapply(names(PR), function(k) setnames(PR[[k]][, .(gene, rank)], "rank", k)))
  sh <- W[Brix <= 100 & NDF <= 100 & NFC <= 100]
  loc <- PR$Brix[gene %in% sh$gene, .(gene, chr, pos)]
  put("nBapNshared", nrow(sh))
  put("nBapNsharedLocus", sum(loc$chr == 4 & loc$pos >= LFROM & loc$pos <= LTO))
  ## the beta-alanine pathway genes' MAGMA p for Brix (other than the ALDH)
  pwf <- Sys.getenv("PATHWAY_FILE", file.path(dirname(BASE), "MAGCAT/inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))
  if (file.exists(pwf)) {
    pw <- fread(pwf, sep = "\t", quote = "")
    bg <- unique(pw[`Pathway-id` == "PWY-3981" & grepl("^SORBI", `Gene-name`), `Gene-name`])
    o  <- PR$Brix[gene %in% setdiff(bg, ALDH), magma_p]
    put("nBapBalaOtherLo", pf(min(o))); put("nBapBalaOtherHi", pf(max(o)))
  }
}

## ---- the locus: leads, LD, genotype classes, allele effects -----------------------
if (file.exists(file.path(LOC, "leads.tsv"))) {
  ld <- fread(file.path(LOC, "leads.tsv"), header = FALSE, col.names = c("trait", "rs", "ps", "p"))
  for (k in names(TR)) { r <- ld[trait == TR[k]]
    if (nrow(r)) { put(sprintf("nBapLeadP%s", k), pf(r$p)); put(sprintf("nBapLeadPos%s", k), sprintf("%.3f", r$ps / 1e6)) } }
  if (file.exists(file.path(LOC, "leads_ld.ld"))) {
    L <- as.matrix(fread(file.path(LOC, "leads_ld.ld"), header = FALSE)); put("nBapLeadRsq", f2(min(L[upper.tri(L)])))
  }
  foc <- ld[trait == TR["Brix"], rs]
  if (!is.null(raw) && length(foc)) {
    gc <- names(raw)[startsWith(names(raw), paste0(foc, "_"))]; d <- raw[[gc]]
    a1 <- substring(gc, nchar(foc) + 2)
    e <- fread(file.path(LOC, "lead_effects.tsv"))[rs == foc]
    up <- with(e[trait == TR["Brix"]], if (beta > 0) allele1 else allele0)
    nup <- if (a1 == up) d else 2 - d
    put("nBapUp", up); put("nBapNUpHom", sum(nup == 2, na.rm = TRUE)); put("nBapNOtherHom", sum(nup == 0, na.rm = TRUE))
    put("nBapNHet", sum(nup == 1, na.rm = TRUE)); put("nBapFreqUp", f2(mean(nup, na.rm = TRUE) / 2))
    for (k in names(TR)) { r <- e[trait == TR[k]]; if (!nrow(r)) next
      s <- if (r$allele1 == up) 1 else -1; sdv <- sd(ph[[TR[k]]], na.rm = TRUE)
      put(sprintf("nBapEff%s", k), f2(s * r$beta / sdv)); put(sprintf("nBapEffP%s", k), pf(r$p_wald)) }
  }
}
sf <- file.path(LOC, "43_chr4_locus_summary.csv")
if (file.exists(sf)) { s <- fread(sf)
  put("nBapNboth", s$n_lines); put("nBapDiffBrix", sprintf("%.1f", 2 * s$b_x)); put("nBapDiffNDF", sprintf("%.1f", 2 * abs(s$b_y)))
  put("nBapAdjNDF", sprintf("%.1f", 2 * abs(s$b_y_adj_x))); put("nBapAdjP", pf(s$p_y_adj_x))
  put("nBapFracVia", sprintf("%.0f", 100 * s$frac_y_via_x)); put("nBapZcor", f2(s$cor_z)); put("nBapNwin", big(s$n_window_snps)) }

## ---- write the macros ---------------------------------------------------------------
ALL <- c("nBapNlines", paste0("nBapN", names(TR)), "nBapRBrixNDF", "nBapRBrixNFC", "nBapRNDFNFC", "nBapRBrixADF",
         "nBapRBrixLig", paste0("nBapGC", names(TR)), paste0("nBapSig", names(TR)), paste0("nBapComp", names(TR)),
         paste0("nBapSigList", names(TR)), "nBapNpw", paste0("nBapBala", names(TR)), paste0("nBapBalaComp", names(TR)),
         "nBapBalaN", paste0("nBapTrxPw", names(TR)), paste0("nBapTrxPwComp", names(TR)), "nBapNgenes",
         paste0("nBapNthree", names(TR)),
         as.vector(outer(c("nBapAldh", "nBapAldhTl", "nBapAldhGwas", "nBapAldhMagmaP",
                           "nBapTrx", "nBapTrxTl", "nBapTrxGwas", "nBapTrxMagmaP"), names(TR), paste0)),
         paste0("nBapTopLoc", names(TR)), paste0("nBapTopP", names(TR)), paste0("nBapLocP", names(TR)),
         "nBapNshared", "nBapNsharedLocus", "nBapBalaOtherLo", "nBapBalaOtherHi",
         paste0("nBapLeadP", names(TR)), paste0("nBapLeadPos", names(TR)), "nBapLeadRsq",
         "nBapUp", "nBapNUpHom", "nBapNOtherHom", "nBapNHet", "nBapFreqUp",
         paste0("nBapEff", names(TR)), paste0("nBapEffP", names(TR)),
         "nBapNboth", "nBapDiffBrix", "nBapDiffNDF", "nBapAdjNDF", "nBapAdjP", "nBapFracVia", "nBapZcor", "nBapNwin")
ALL <- unique(ALL)
writeLines(c("%% Placeholders for the macros written by 44_bap_chr4_numbers.R",
             sprintf("\\providecommand{\\%s}{\\textbf{??}}", ALL)), file.path(TEX, "bap_chr4_stub.tex"))
have <- intersect(ALL, names(M))
writeLines(c("%% written by 44_bap_chr4_numbers.R -- do not edit",
             sprintf("\\renewcommand{\\%s}{%s}", have, unlist(M[have]))), file.path(TEX, "bap_chr4.tex"))
cat(sprintf("%d of %d macros set; still ??: %s\n", length(have), length(ALL),
            paste(setdiff(ALL, have), collapse = ", ")))

## ---- table: genes of the locus --------------------------------------------------------
if (length(PR) == 3) {
  ANN <- c(SORBI_3004G300600 = "putative receptor-like kinase",
           SORBI_3004G300800 = "\\textit{SbALDH3E1}, ALDH3 fatty aldehyde dehydrogenase",
           SORBI_3004G300900 = "NADPH-thioredoxin reductase",
           SORBI_3004G301100 = "GATA transcription factor")
  B <- PR$Brix[as.character(chr) == "4" & pos >= LFROM & pos <= LTO, .(gene, pos, best_pathway_name)]
  for (k in names(TR)) B <- merge(B, PR[[k]][, .(gene, r = rank, g = gwas_rank)][, setnames(.SD, c("r", "g"), paste0(c("r", "g"), k))], by = "gene")
  B <- B[pmin(rBrix, rNDF, rNFC) <= 100][order(pos)]
  cell <- function(r, g) sprintf("%s (%s)", rk(r), rk(g))
  rows <- B[, sprintf("Sobic.%s & %.3f & %s & %s & %s & %s & %s \\\\",
                      sub("SORBI_3", "", gene), pos / 1e6,
                      ifelse(gene %in% names(ANN), ANN[gene], "--"),
                      ifelse(is.na(best_pathway_name) | best_pathway_name == "", "--", tex_esc(pw_name(best_pathway_name))),
                      mapply(cell, rBrix, gBrix), mapply(cell, rNDF, gNDF), mapply(cell, rNFC, gNFC))]
  out <- c("\\begin{table*}[!t]", "\\centering",
           "\\caption{Genes of the chromosome~4 locus shared by stem sugar and fibre in the sorghum BAP (GEMMA with three principal components, MAGMA, CATFISH). Genes with a gene-score rank of $100$ or better for at least one trait, in order of position (gene start, BTx623 v3). Each trait column gives the CATFISH rank among all genes, and in parentheses the rank on the GWAS layer alone (best variant $p$-value; genes that share their best variant share this rank, given as the mid-rank rounded up). Pathway: the gene's best-ranked SorghumBicolorCyc pathway, which is what moves a gene above neighbours that share its GWAS evidence; -- no pathway. Annotation from SorghumBicolorCyc, the sorghum ALDH superfamily \\citep{Islam2022SbALDH} and PlantTFDB/PRGDB; -- not examined.}",
           "\\label{tab:bap_chr4}", "\\small",
           "\\begin{tabular}{lrp{4.4cm}p{3.6cm}rrr}", "\\hline",
           "\\textbf{Gene} & \\textbf{Mb} & \\textbf{Annotation} & \\textbf{Pathway} & \\textbf{Brix} & \\textbf{NDF} & \\textbf{NFC} \\\\",
           "\\hline", rows, "\\hline", "\\end{tabular}", "\\end{table*}")
  writeLines(out, file.path(TEX, "tab_bap_chr4.tex"))
  cat("wrote", file.path(TEX, "tab_bap_chr4.tex"), sprintf("(%d genes)\n", nrow(B)))
}
cat("wrote", file.path(TEX, "bap_chr4.tex"), "\n")
