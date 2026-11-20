#!/usr/bin/env Rscript
## =============================================================================
## 38_competitive_real.R -- CATFISH's self-contained and competitive tests side
## by side in the real analyses: the ten GLGC lipid sets and sorghum dry
## biomass (Section 3.11, Table tab:competitive_real).
##
## Reads  catfish2/results/gc_real_comp/<set>_..._5comp_westfall_young.csv, the
##        COMP=1 runs of scripts/hpc/gc_real/wy_gc_real.R (centred gene p-values,
##        B = 10,000 for both tests); their self-contained columns are those
##        of the paper's runs (checked in 36_prioritize_comp.R)
## Writes catfish2/tex/auto/tab_competitive_real.tex
##        catfish2/results/38_competitive_real.csv
## Calls: Westfall-Young family-wise p <= 0.05 (p_fwer_std, p_fwer_comp); both:
##        max of the two <= 0.05. Size: Spearman correlation between pathway
##        size and evidence (-z), as in tab_magma_compare. Agreement: Spearman
##        correlation between the two standardised statistics.
## Usage (from catfish2/scripts/final):  Rscript 38_competitive_real.R
## =============================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
CMP  <- Sys.getenv("COMP_DIR", file.path(BASE, "catfish2/results/gc_real_comp"))
AUT  <- Sys.getenv("TEX_AUTO", file.path(BASE, "catfish2/tex/auto"))
RES  <- Sys.getenv("RES_OUT",  file.path(BASE, "catfish2/results"))
LIPO <- "REACTOME_PLASMA_LIPOPROTEIN_ASSEMBLY_REMODELING_AND_CLEARANCE"

SETS <- data.table(
  set   = c(c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_")),
            "Dry_tons_per_acre"),
  trait = c(rep(c("LDL-C", "HDL-C"), 5), "Dry biomass"),
  group = c(rep(c("EUR", "AFR", "EAS", "SAS", "HIS"), each = 2), "Sorghum"))
SETS[, file := ifelse(set == "Dry_tons_per_acre",
                      sprintf("%s_P_SNPWISE_MEAN_5comp_westfall_young.csv", set),
                      sprintf("%s_reactome_5comp_westfall_young.csv", set))]

sp <- function(a, b) suppressWarnings(stats::cor(a, b, method = "spearman", use = "complete.obs"))
rows <- lapply(seq_len(nrow(SETS)), function(i) {
  f <- file.path(CMP, SETS$file[i]); if (!file.exists(f)) stop("missing: ", f, call. = FALSE)
  w <- fread(f)
  rk <- function(z) frank(z, ties.method = "min")
  lip <- w$pathway_id == LIPO
  data.table(SETS[i, .(set, trait, group)], K = nrow(w),
             sig_self = sum(w$p_fwer_std <= 0.05), sig_comp = sum(w$p_fwer_comp <= 0.05),
             sig_both = sum(pmax(w$p_fwer_std, w$p_fwer_comp) <= 0.05),
             lipo_rank_self = if (any(lip)) rk(w$z_std)[lip] else NA_integer_,
             lipo_rank_comp = if (any(lip)) rk(w$z_std_comp)[lip] else NA_integer_,
             top_comp = w$pathway_id[which.min(w$z_std_comp)],
             agree = sp(w$z_std, w$z_std_comp),
             size_self = sp(w$n_genes, -w$z_std), size_comp = sp(w$n_genes, -w$z_std_comp))
})
D <- rbindlist(rows)
fwrite(D, file.path(RES, "38_competitive_real.csv"))

big <- function(x) formatC(x, format = "d", big.mark = "{,}")
f2  <- function(x) sub("^\\$-0\\.00", "$0.00", sprintf("$%.2f$", x))
na  <- function(x) ifelse(is.na(x), "--", as.character(x))
L <- c("\\begin{table*}[!t]", "\\centering",
  paste0("\\caption{CATFISH's self-contained and competitive tests on the same centred gene ",
         "statistics in the real analyses ($1{,}390$ Reactome pathways for the lipid sets, ",
         "$", D[set == "Dry_tons_per_acre", K], "$ SorghumBicolorCyc pathways for sorghum; ",
         "$B=10{,}000$ for each test). Significant: Westfall--Young family-wise $p\\le0.05$; ",
         "both: significant in both tests. Lipoprotein rank: rank of plasma lipoprotein ",
         "assembly, remodeling and clearance. Agreement: Spearman correlation between the two ",
         "standardised statistics. Size: Spearman correlation between pathway size and ",
         "pathway evidence ($-z_{\\mathrm{std}}$).}"),
  "\\label{tab:competitive_real}", "\\small",
  "\\begin{tabular}{llrrrccrcrr}", "\\hline",
  paste0(" & & \\multicolumn{3}{c}{\\textbf{Significant}} & \\multicolumn{2}{c}{\\textbf{Lipoprotein rank}} ",
         "& & & \\multicolumn{2}{c}{\\textbf{Size}} \\\\"),
  paste0("\\textbf{Trait} & \\textbf{Group} & Self & Comp. & Both & Self & Comp. & ",
         "\\textbf{Agreement} & & Self & Comp. \\\\"), "\\hline")
for (i in seq_len(nrow(D))) {
  r <- D[i]
  if (r$set == "Dry_tons_per_acre") L <- c(L, "\\hline")
  L <- c(L, sprintf("%s & %s & %s & %s & %s & %s & %s & %s & & %s & %s \\\\",
    r$trait, r$group, big(r$sig_self), big(r$sig_comp), big(r$sig_both),
    na(r$lipo_rank_self), na(r$lipo_rank_comp), f2(r$agree), f2(r$size_self), f2(r$size_comp)))
}
L <- c(L, "\\hline", "\\end{tabular}", "\\end{table*}")
writeLines(L, file.path(AUT, "tab_competitive_real.tex"))
print(D[, .(set, sig_self, sig_comp, sig_both, lipo_rank_self, lipo_rank_comp,
            agree = round(agree, 2), size_self = round(size_self, 2),
            size_comp = round(size_comp, 2), top_comp = substr(top_comp, 10, 45))], row.names = FALSE)
cat("wrote", file.path(AUT, "tab_competitive_real.tex"), "\n")
