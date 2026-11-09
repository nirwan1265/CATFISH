#!/usr/bin/env Rscript
## =============================================================================
## 34_compare_magma_lipids.R -- CATFISH against MAGMA's own gene-set analysis
## (competitive, and --model self-contained) on the ten GLGC lipid sets, over
## the 1,390 Reactome pathways CATFISH tested (10-200 genes).
## Inputs:  results/human_LDL/catfish/<ta>_reactome_5comp_westfall_young.csv (CATFISH, GC)
##          catfish2/results/magma_gsa/<ta>.gsa.out, <ta>.gsa.self.out (MAGMA v1.10,
##          merged .genes.raw, the same Reactome .sets.annot)
## Writes   catfish2/results/magma_gsa/34_compare_magma_lipids.csv
##          catfish2/tex/auto/tab_magma_compare.tex
##   cat_fwer   CATFISH Westfall-Young family-wise p <= 0.05
##   mc_bonf / mc_bh   MAGMA competitive: Bonferroni p <= 0.05/K, BH q <= 0.05
##   ms_bonf    MAGMA self-contained, Bonferroni
##   lp_*       rank of plasma lipoprotein assembly, remodeling and clearance
##   fam10_*    plasma-lipoprotein subtree sets among the top ten
##   rho_*      Spearman correlation of pathway evidence with CATFISH
##   size_*     Spearman correlation of pathway evidence with pathway size
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
MG <- Sys.getenv("MAGMA_GSA_DIR", file.path(BASE, "catfish2/results/magma_gsa"))
CF <- Sys.getenv("LIPID_WY_DIR", file.path(BASE, "results/human_LDL/catfish"))
TEX <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
TA <- c("ldl_eur","hdl_eur","ldl_afr","hdl_afr","ldl_eas","hdl_eas","ldl_sas","hdl_sas","ldl_his","hdl_his")
LP <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"
RX <- "^REACTOME_(PLASMA_LIPOPROTEIN|LDL_|HDL_|VLDL_|VLDLR_|CHYLOMICRON_)"
rd <- function(f) { ln <- readLines(f); ln <- ln[!startsWith(ln, "#")]; d <- fread(text = paste(ln, collapse = "\n"), fill = TRUE)
  id <- if ("FULL_NAME" %in% names(d)) ifelse(is.na(d$FULL_NAME) | d$FULL_NAME == "", d[[1]], d$FULL_NAME) else d[[1]]
  data.table(pathway_id = id, ngenes = d$NGENES, p = d$P) }
rows <- list(); tops <- list()
for (ta in TA) {
  cf <- fread(file.path(CF, paste0(ta, "_reactome_5comp_westfall_young.csv")))
  K <- nrow(cf)
  mc <- rd(file.path(MG, paste0(ta, ".gsa.out")))[pathway_id %in% cf$pathway_id]
  ms <- rd(file.path(MG, paste0(ta, ".gsa.self.out")))[pathway_id %in% cf$pathway_id]
  setorder(cf, z_std); setorder(mc, p); setorder(ms, p)
  m <- merge(cf[, .(pathway_id, n_genes, z_std, p_fwer_std)], mc[, .(pathway_id, p_comp = p)], by = "pathway_id")
  m <- merge(m, ms[, .(pathway_id, p_self = p)], by = "pathway_id")
  sp <- function(a, b) round(cor(a, b, method = "spearman", use = "complete.obs"), 2)
  rows[[ta]] <- data.table(set = ta, K = K, matched = nrow(m),
    cat_fwer = sum(cf$p_fwer_std <= 0.05),
    mc_bonf = sum(mc$p <= 0.05 / K), mc_bh = sum(p.adjust(mc$p, "BH") <= 0.05),
    ms_bonf = sum(ms$p <= 0.05 / K),
    lp_cat = match(LP, cf$pathway_id), lp_mc = match(LP, mc$pathway_id), lp_ms = match(LP, ms$pathway_id),
    lp_mc_p = signif(mc$p[match(LP, mc$pathway_id)], 2),
    fam10_cat = sum(grepl(RX, cf$pathway_id[1:10])), fam10_mc = sum(grepl(RX, mc$pathway_id[1:10])),
    rho_cat_mc = sp(-m$z_std, -log10(m$p_comp)), rho_cat_ms = sp(-m$z_std, -log10(m$p_self)),
    size_cat = sp(m$n_genes, -m$z_std), size_mc = sp(m$n_genes, -log10(m$p_comp)), size_ms = sp(m$n_genes, -log10(m$p_self)),
    top20_overlap = length(intersect(cf$pathway_id[1:20], mc$pathway_id[1:20])))
  if (ta %in% c("ldl_eur","hdl_eur")) tops[[ta]] <- data.table(set = ta, rank = 1:10,
    catfish = substr(sub("^REACTOME_", "", cf$pathway_id[1:10]), 1, 45),
    magma_comp = substr(sub("^REACTOME_", "", mc$pathway_id[1:10]), 1, 45), p_comp = signif(mc$p[1:10], 2))
}
r <- rbindlist(rows); print(r[, 1:12], row.names = FALSE); cat("\n"); print(r[, c(1, 13:18)], row.names = FALSE)
cat("\n"); print(rbindlist(tops), row.names = FALSE)
fwrite(r, file.path(MG, "34_compare_magma_lipids.csv"))

## ---- LaTeX table ---------------------------------------------------------------
lab <- c(eur = "EUR", afr = "AFR", eas = "EAS", sas = "SAS", his = "HIS")
num <- function(x) formatC(x, format = "d", big.mark = "{,}")
f2  <- function(x) { x <- round(x, 2); x[x == 0] <- 0; formatC(x, format = "f", digits = 2) }
L <- c("\\begin{table*}[!t]", "\\centering",
  "\\caption{CATFISH and MAGMA's gene-set analysis on the same gene-level statistics and the same $1{,}390$ Reactome pathways (10--200 genes), in the ten GLGC lipid sets. Significant: CATFISH Westfall--Young family-wise $p\\le0.05$; MAGMA Bonferroni $p\\le0.05/1{,}390$. Lipoprotein rank: rank of plasma lipoprotein assembly, remodeling and clearance. Size: Spearman correlation between pathway size and pathway evidence (CATFISH $-z_{\\mathrm{std}}$, MAGMA $-\\log_{10}p$). CATFISH uses centred gene $p$-values (Section~2.7); MAGMA reads the uncentred gene statistics.}",
  "\\label{tab:magma_compare}", "\\small",
  "\\begin{tabular}{llrrrccccrrr}", "\\hline",
  " & & \\multicolumn{3}{c}{\\textbf{Significant}} & \\multicolumn{3}{c}{\\textbf{Lipoprotein rank}} & & \\multicolumn{3}{c}{\\textbf{Size}} \\\\",
  "\\textbf{Trait} & \\textbf{Group} & CATFISH & \\makecell[r]{MAGMA\\\\self} & \\makecell[r]{MAGMA\\\\comp.} & CATFISH & \\makecell[c]{MAGMA\\\\self} & \\makecell[c]{MAGMA\\\\comp.} & & CATFISH & \\makecell[r]{MAGMA\\\\self} & \\makecell[r]{MAGMA\\\\comp.} \\\\", "\\hline")
for (i in seq_len(nrow(r))) {
  x <- r[i]; tr <- toupper(sub("_.*", "", x$set)); an <- lab[[sub(".*_", "", x$set)]]
  L <- c(L, sprintf("%s & %s & %s & %s & %s & %d & %d & %d & & $%s$ & $%s$ & $%s$ \\\\",
    if (tr == "LDL") "LDL-C" else "HDL-C", an, num(x$cat_fwer), num(x$ms_bonf), num(x$mc_bonf),
    x$lp_cat, x$lp_ms, x$lp_mc, f2(x$size_cat), f2(x$size_ms), f2(x$size_mc)))
}
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(TEX, "tab_magma_compare.tex"))
cat("wrote", file.path(TEX, "tab_magma_compare.tex"), "\n")
