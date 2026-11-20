#!/usr/bin/env Rscript
## =============================================================================
## 36_prioritize_comp.R -- the gene prioritisation of Section 2.5 with the
## pathway layer taken from the competitive test, or from the joint call,
## instead of the self-contained test. Everything else (GWAS and MAGMA layers,
## weights, hub term, pathway sets) is unchanged; 06 and 25 are run as they are,
## only with a different pathway statistic and output folder.
##
## Pathway statistic (smaller = stronger), from the COMP=1 Westfall-Young runs
## (scripts/hpc/gc_real/wy_gc_real.R, out_comp/), which carry both tests:
##   self   z_std        the paper's runs; read from their existing outputs
##   comp   z_std_comp   competitive (circular genomic permutation)
##   joint  z_std_joint  = max(z_std, z_std_comp), the weaker of the two
##                        standardised statistics: a pathway ranks high only
##                        if it is strong on both tests (as the joint call)
##
## Inputs
##   catfish2/results/gc_real_comp/<trait>_<anc>_reactome_5comp_westfall_young.csv
##   catfish2/results/gc_real_comp/<TRAIT>_P_SNPWISE_MEAN_5comp_westfall_young.csv
##   (rsync of gc_real/out_comp/ from the cluster). The self-contained z_std in
##   these files is checked against the file the paper used; the script stops
##   if they differ.
## Writes, inside catfish2/results/prioritize_comp/:
##   wy/        the competitive WY files with the z_std_joint column added
##   comp/      06_candidate_genes_<tag>.csv, sorghum/25_<TRAIT>_gene_prioritize.csv
##   joint/     the same with the joint statistic
## Runs that already have their output are skipped (FORCE=1 reruns them).
##
## Usage (from catfish2/scripts):  Rscript 36_prioritize_comp.R
##   then                          Rscript 37_prioritize_comp_recovery.R
## Env: CATFISH_BASE, COMP_DIR, VARIANTS (default "comp,joint"),
##      SETS_RUN (default all 12, e.g. "ldl_eur,hdl_eur"), FORCE
## =============================================================================
suppressPackageStartupMessages(library(data.table))

BASE     <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
COMP_DIR <- Sys.getenv("COMP_DIR", file.path(BASE, "catfish2/results/gc_real_comp"))
OUT      <- file.path(BASE, "catfish2/results/prioritize_comp")
VARIANTS <- strsplit(Sys.getenv("VARIANTS", "comp,joint"), ",")[[1]]
FORCE    <- Sys.getenv("FORCE", "0") == "1"
STAT     <- c(comp = "z_std_comp", joint = "z_std_joint")

LIPID   <- c(outer(c("ldl", "hdl"), c("eur", "afr", "eas", "sas", "his"), paste, sep = "_"))
SORGHUM <- c("Dry_tons_per_acre", "Lignin_DM_perc")
ALL     <- c(LIPID, SORGHUM)
SETS_RUN <- strsplit(Sys.getenv("SETS_RUN", paste(ALL, collapse = ",")), ",")[[1]]
stopifnot(all(SETS_RUN %in% ALL), all(VARIANTS %in% names(STAT)))

wy_name  <- function(s) if (s %in% SORGHUM) sprintf("%s_P_SNPWISE_MEAN_5comp_westfall_young.csv", s) else
  sprintf("%s_reactome_5comp_westfall_young.csv", s)
paper_wy <- function(s) if (s %in% SORGHUM) file.path(BASE, "catfish2/results/sorghum_bap", wy_name(s)) else
  file.path(BASE, "results/human_LDL/catfish", wy_name(s))

## ---- 1. competitive WY files, checked, with the joint statistic added ---------
dir.create(file.path(OUT, "wy"), recursive = TRUE, showWarnings = FALSE)
miss <- SETS_RUN[!file.exists(file.path(COMP_DIR, vapply(SETS_RUN, wy_name, "")))]
if (length(miss)) stop("no competitive WY result for: ", paste(miss, collapse = ", "),
                       "\n  expected in ", COMP_DIR, call. = FALSE)
for (s in SETS_RUN) {
  w <- fread(file.path(COMP_DIR, wy_name(s)))
  p <- fread(paper_wy(s))
  if (!all(c("z_std", "z_std_comp") %in% names(w)))
    stop(wy_name(s), ": needs z_std and z_std_comp (a COMP=1 run)", call. = FALSE)
  d <- max(abs(w$z_std - p$z_std[match(w$pathway_id, p$pathway_id)]), na.rm = TRUE)
  if (!is.finite(d) || d > 1e-6 || nrow(w) != nrow(p))
    stop(sprintf("%s: self-contained z_std differs from the paper's file (max diff %g)", s, d),
         call. = FALSE)
  w[, z_std_joint := pmax(z_std, z_std_comp)]
  fwrite(w, file.path(OUT, "wy", wy_name(s)))
  cat(sprintf("%-20s %4d pathways; self z_std matches the paper's file\n", s, nrow(w)))
}

## ---- 2. rerun 06 (lipids) and 25 (sorghum) with the other pathway statistic ---
run <- function(script, env) {
  cat(sprintf("\n>>> %s  %s\n", script, paste(env, collapse = " ")))
  st <- system2("Rscript", script, env = env)
  if (st != 0) stop(script, " failed (exit ", st, ") for ", paste(env, collapse = " "), call. = FALSE)
}
for (v in VARIANTS) {
  vdir <- file.path(OUT, v)
  dir.create(file.path(vdir, "sorghum"), recursive = TRUE, showWarnings = FALSE)
  for (s in intersect(SETS_RUN, LIPID)) {
    tr <- sub("_.*", "", s); anc <- sub(".*_", "", s)
    tag <- if (anc == "eur") tr else s
    if (!FORCE && file.exists(file.path(vdir, sprintf("06_candidate_genes_%s.csv", tag)))) {
      cat(sprintf("skip %s %s (done)\n", v, s)); next }
    run("06_candidate_genes.R",
        c("SPECIES=human", paste0("TRAIT=", tr), paste0("ANC=", anc),
          paste0("WY_CSV=", file.path(OUT, "wy", wy_name(s))),
          paste0("PATHWAY_P_COL=", STAT[[v]]),
          paste0("CANDIDATE_OUT_DIR=", vdir)))
    ## 06 sets a layer to zero when it cannot find its input; the comparison is
    ## only fair if every layer covers the same genes as in the paper's run
    a <- fread(file.path(BASE, "catfish2/results", sprintf("06_candidate_genes_%s.csv", tag)))
    b <- fread(file.path(vdir, sprintf("06_candidate_genes_%s.csv", tag)))
    for (cc in c("rank_gwas", "rank_magma", "rank_pathway"))
      if (sum(!is.na(a[[cc]])) != sum(!is.na(b[[cc]])))
        stop(sprintf("%s %s: %s covers %d genes here and %d in the paper's run (set GWAS_FILE?)",
                     v, s, cc, sum(!is.na(b[[cc]])), sum(!is.na(a[[cc]]))), call. = FALSE)
  }
  sg <- intersect(SETS_RUN, SORGHUM)
  sg <- sg[FORCE | !file.exists(file.path(vdir, "sorghum", sprintf("25_%s_gene_prioritize.csv", sg)))]
  if (length(sg))
    run("25_sorghum_gene_prioritize.R",
        c(paste0("TRAITS=", paste(sg, collapse = ",")),
          paste0("WY_DIR=", file.path(OUT, "wy")),
          paste0("PATHWAY_STAT_COL=", STAT[[v]]),
          paste0("OUT_DIR=", file.path(vdir, "sorghum"))))
}
cat("\ndone; next: Rscript 37_prioritize_comp_recovery.R\n")
