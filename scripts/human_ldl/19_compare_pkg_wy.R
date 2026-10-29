## =============================================================================
## 19_compare_pkg_wy.R
##
## Checks that the CATFISH package's Westfall-Young engine reproduces
## 18_westfall_young_trait.R. Same inputs, same gene filter, same pathway
## parsing, same size window, same seed, same B and TFisher setting, so the
## two must agree to the last digit.
##
##   route "wy"    -> catfish_wy_pathways()            (stand-alone function)
##   route "omni2" -> catfish_omni2_pathways(perm_mode = "wy")
##   route "both"  -> both (runs the null twice)
##
## USAGE (B and use_tfisher MUST match the 18 run being compared):
##   TRAIT=ldl Rscript 19_compare_pkg_wy.R 10000 TRUE both
##   TRAIT=hdl Rscript 19_compare_pkg_wy.R 10000 TRUE wy
##
## Reads   results/human_LDL/catfish/<trait>_eur_reactome_<k>comp_westfall_young.csv
## Writes  results/human_LDL/catfish/<trait>_eur_reactome_<k>comp_wy_pkg.csv
##         results/human_LDL/catfish/<trait>_eur_reactome_<k>comp_omni2_wy_pkg.csv
## =============================================================================
suppressPackageStartupMessages({ library(CATFISH) })

args <- commandArgs(trailingOnly = TRUE)
B           <- if (length(args) >= 1) as.integer(args[[1]]) else 1000L
USE_TFISHER <- if (length(args) >= 2) as.logical(args[[2]]) else FALSE
ROUTE       <- if (length(args) >= 3) args[[3]] else "both"
if (!ROUTE %in% c("wy", "omni2", "both")) stop("route must be wy, omni2 or both", call. = FALSE)

## ---- identical to the human branch of 18_westfall_young_trait.R
SEED      <- 20260706
MIN_GENES <- 10L
MAX_GENES <- 200L
BASE   <- path.expand("~/Documents/Github/MAGCAT")
DATA_D <- path.expand("~/Documents/Github/MAGCAT_data")
TRAIT  <- tolower(Sys.getenv("TRAIT", "ldl"))
MDIR <- if (TRAIT == "ldl") file.path(BASE, "results/human_LDL/magma") else
        file.path(BASE, "results/human_LDL", paste0("magma_", TRAIT))
GENES_FILE <- file.path(MDIR, sprintf("%s_eur.genes.out", TRAIT))
CORPAIRS   <- file.path(BASE, "results/human_LDL/catfish", sprintf("%s_eur.cor_pairs.txt", TRAIT))
SETS       <- file.path(DATA_D, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")
OUT_D      <- file.path(BASE, "results/human_LDL/catfish")
TAG        <- sprintf("%s_eur_reactome_%dcomp", TRAIT, if (isTRUE(USE_TFISHER)) 5L else 4L)

F18 <- file.path(OUT_D, paste0(TAG, "_westfall_young.csv"))
if (!file.exists(F18)) stop("18 output not found: ", F18, call. = FALSE)
r18 <- read.csv(F18, stringsAsFactors = FALSE)
cat(sprintf("trait %s | B = %d | TFisher = %s | route = %s\n", TRAIT, B, USE_TFISHER, ROUTE))
cat(sprintf("18 output: %s (%d pathways)\n", basename(F18), nrow(r18)))

## the 18 output does not record B: its smallest attainable WY p is 1/(B+1)
B18 <- round(1 / min(r18$p_fwer_wy, na.rm = TRUE)) - 1
if (B18 > B) stop(sprintf("the 18 file was run with B >= %d but you asked for B = %d", B18, B),
                  call. = FALSE)

## ---- inputs, read and filtered exactly as 18 does
g <- read.table(GENES_FILE, header = TRUE, comment.char = "#", stringsAsFactors = FALSE)
g$GENE <- as.character(g$GENE)
g <- g[is.finite(g$P) & is.finite(g$ZSTAT), ]
lines <- readLines(SETS, warn = FALSE)
pw <- list()
for (ln in lines) {
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (length(f) < 2) next
  gs <- unique(f[-1]); pw[[f[1]]] <- gs[nzchar(gs)]
}
cp <- if (requireNamespace("data.table", quietly = TRUE)) {
  as.data.frame(data.table::fread(CORPAIRS, showProgress = FALSE))
} else {
  utils::read.delim(CORPAIRS, header = TRUE, stringsAsFactors = FALSE,
                    colClasses = c("character", "character", "numeric"))
}

## ---- comparison
compare <- function(pk, label) {
  cat(sprintf("\n==== %s vs 18 ====\n", label))
  only18 <- setdiff(r18$pathway_id, pk$pathway_id)
  onlypk <- setdiff(pk$pathway_id, r18$pathway_id)
  cat(sprintf("pathways: 18 = %d, package = %d, only in 18 = %d, only in package = %d\n",
              nrow(r18), nrow(pk), length(only18), length(onlypk)))
  cols <- c("z_std", "p_fwer_std", "n_genes", "log10_p_omnibus", "p_omnibus",
            "p_fwer_wy", "p_bonferroni")
  # 18's values went through write.csv / read.csv: send the package's the same way
  tmp <- tempfile(fileext = ".csv")
  write.csv(pk[, c("pathway_id", cols)], tmp, row.names = FALSE)
  pk <- read.csv(tmp, stringsAsFactors = FALSE)
  m <- match(r18$pathway_id, pk$pathway_id)
  ok <- !is.na(m)
  out <- data.frame(column = cols, identical = NA_integer_, n = sum(ok),
                    max_abs_diff = NA_real_, max_rel_diff = NA_real_)
  for (j in seq_along(cols)) {
    a <- r18[[cols[j]]][ok]; b <- pk[[cols[j]]][m[ok]]
    same <- (is.na(a) & is.na(b)) | (!is.na(a) & !is.na(b) & a == b)
    d <- abs(a - b)
    out$identical[j]    <- sum(same)
    out$max_abs_diff[j] <- if (any(is.finite(d))) max(d[is.finite(d)]) else 0
    rel <- d / pmax(abs(a), 1e-300)
    out$max_rel_diff[j] <- if (any(is.finite(rel))) max(rel[is.finite(rel)]) else 0
  }
  print(out, row.names = FALSE)
  n_sig18 <- sum(r18$p_fwer_std <= 0.05, na.rm = TRUE)
  n_sigpk <- sum(pk$p_fwer_std <= 0.05, na.rm = TRUE)
  cat(sprintf("significant at FWER 5%% (WY standardised): 18 = %d, package = %d\n",
              n_sig18, n_sigpk))
  agree <- all(out$identical == out$n) && !length(only18) && !length(onlypk)
  cat(if (agree) "RESULT: IDENTICAL\n" else
        "RESULT: DIFFERENT -- see max_abs_diff / max_rel_diff above\n")
  invisible(out)
}

if (ROUTE %in% c("wy", "both")) {
  t0 <- Sys.time()
  wy <- catfish_wy_pathways(
    g, pathways = pw, gene_col = "GENE", p_col = "P", chr_col = "CHR",
    z_col = "ZSTAT", cor_pairs = cp, B = B, use_tfisher = USE_TFISHER,
    min_genes = MIN_GENES, max_genes = MAX_GENES, seed = SEED)
  cat(sprintf("catfish_wy_pathways: %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  write.csv(wy, file.path(OUT_D, paste0(TAG, "_wy_pkg.csv")), row.names = FALSE)
  compare(wy, "catfish_wy_pathways()")
}

if (ROUTE %in% c("omni2", "both")) {
  t0 <- Sys.time()
  om <- catfish_omni2_pathways(
    g, pathways = pw, gene_col = "GENE", p_raw_col = "P", z_col = "ZSTAT",
    chr_col = "CHR", min_genes = MIN_GENES, max_genes = MAX_GENES,
    use_tfisher = USE_TFISHER, perm_mode = "wy", B_perm = B,
    magma_cor_pairs = cp, seed = SEED)
  cat(sprintf("catfish_omni2_pathways(perm_mode = 'wy'): %.1f min\n",
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  write.csv(om, file.path(OUT_D, paste0(TAG, "_omni2_wy_pkg.csv")), row.names = FALSE)
  om_std <- data.frame(
    pathway_id      = om$pathway_id,
    z_std           = om$omni_z_std,
    p_fwer_std      = om$omni_p_fwer_wy_std,
    n_genes         = om$n_genes,
    log10_p_omnibus = om$omni_log10_p_wy_stat,
    p_omnibus       = om$omni_p_wy_stat,
    p_fwer_wy       = om$omni_p_fwer_wy,
    p_bonferroni    = pmin(om$omni_p_wy_stat * nrow(om), 1),
    stringsAsFactors = FALSE)
  compare(om_std, "catfish_omni2_pathways(perm_mode = 'wy')")
}
