#!/usr/bin/env Rscript
## =============================================================================
## 06_candidate_genes.R -- three-layer candidate-gene prioritisation.
##
## Runs on either real dataset:  SPECIES=maize (default) or SPECIES=human.
## Same formula, same weights; only the input paths and ID space differ.
##
## Score(g) = S_GWAS(g) + S_MAGMA(g) + 0.5 * S_PATH(g) + 0.3 * B_hub(g)
##   S(.)   = -log10(rank / N)                    rank 1 = strongest
##   B_hub  = log10(n_pathways(g) + 1)
##
## Pathway evidence is ranked on `z_std`, the Westfall-Young standardised
## statistic (each pathway's observed omnibus against its own null). It gives
## the same ordering as the headline family-wise p-value p_fwer_std, but that
## p-value is floored at 1/(B+1) and ties the whole top of the list, while
## z_std has no floor. (Until 2026-10-29 the default was `log10_p_omnibus`,
## the raw omnibus, whose ordering is dominated by pathway size; set
## PATHWAY_P_COL=log10_p_omnibus to reproduce that.)
##
## Every input path can be overridden with an environment variable so the
## script runs unchanged on another machine.
##
## Usage:  Rscript 06_maize_candidate_genes.R
## =============================================================================

suppressPackageStartupMessages({ library(data.table) })

SPECIES <- tolower(Sys.getenv("SPECIES", "maize"))
if (!SPECIES %in% c("maize", "human"))
  stop("SPECIES must be 'maize' or 'human'", call. = FALSE)
## For human, TRAIT selects which lipid GWAS: ldl (default), hdl, logtg, tc.
## All share the reference panel, the gene annotation and the Reactome sets;
## only the summary statistics and the WY result differ.
TRAIT <- tolower(Sys.getenv("TRAIT", if (SPECIES == "human") "ldl" else SPECIES))
TAG   <- if (SPECIES == "human") TRAIT else SPECIES
BASE <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
DATA_D <- Sys.getenv("CATFISH_DATA", path.expand("~/Documents/Github/MAGCAT_data"))
RES  <- Sys.getenv("CATFISH_RESULTS",
                   path.expand("~/Documents/Research/results/MAGMA/MAGCAT"))
MD   <- file.path(RES, "magma_multi_snp_wise_genes_by_chr_N_maize")

if (SPECIES == "maize") {
  GENES_OUT  <- file.path(MD, "N_maize_MLM_ALLCHR.multi_snp_wise.genes.out")
  WY_CSV     <- file.path(BASE, "catfish2/results/maize_N_MLM_corncyc_westfall_young.csv")
  SETS       <- file.path(BASE, "catfish2/data/maize_pathway.sets.annot.WIDE.txt")
  GENE_LOC   <- file.path(BASE, "catfish2/data/maize.genes.loc")
  PW_NAMES_F <- file.path(BASE, "catfish2/data/corncyc_pathways.20230103.maizeIDs")
  SYMBOL_COL <- NA_integer_          # maize gene IDs are already readable
  GENES_ANNOT <- ""                  # set MAGMA_GENES_ANNOT if you have one
  TARGETS    <- c("Zm00001eb220040") # aspartate kinase named in the prior paper
} else {
  ## LDL sits in magma/, the other traits in magma_<trait>/ beside it; a
  ## couple of other layouts are tolerated so this does not break again
  GENES_OUT <- ""
  for (cand in c(file.path(BASE, "results/human_LDL/magma",
                           sprintf("%s_eur.genes.out", TRAIT)),
                 file.path(BASE, sprintf("results/human_LDL/magma_%s", TRAIT),
                           sprintf("%s_eur.genes.out", TRAIT)),
                 file.path(BASE, sprintf("results/lipids/%s/magma", TRAIT),
                           sprintf("%s_eur.genes.out", TRAIT)))) {
    if (file.exists(cand)) { GENES_OUT <- cand; break }
  }
  if (!nzchar(GENES_OUT))
    stop(sprintf("no %s_eur.genes.out found under results/human_LDL/", TRAIT),
         call. = FALSE)
  ## 5-component results at B = 10000; the superseded B = 200 LDL file has
  ## been deleted, and the 4-component HDL file is kept only as the ablation
  ## control for Section 3.9.3 and is not used here
  WY_CSV     <- file.path(BASE, "results/human_LDL/catfish",
                          sprintf("%s_eur_reactome_5comp_westfall_young.csv", TRAIT))
  SETS       <- file.path(DATA_D, "c2.cp.reactome.v2026.2.Hs.entrez.sets.annot")
  GENE_LOC   <- file.path(DATA_D, "NCBI37.3.gene.loc")
  PW_NAMES_F <- ""                   # Reactome set names are already readable
  SYMBOL_COL <- 6L                   # NCBI37.3.gene.loc: entrez chr start stop strand symbol
  ## the annotation is trait-independent: it maps genes to reference-panel
  ## variants and never sees the GWAS, so one file serves every trait
  GENES_ANNOT <- file.path(DATA_D, "g1000_eur.genes.annot")
  if (!file.exists(GENES_ANNOT))
    GENES_ANNOT <- file.path(DATA_D, "ldl_eur.genes.annot")
  ## canonical monogenic / established LDL-C genes, as a positive control
  TARGETS <- switch(TRAIT,
    ldl   = c("LDLR","PCSK9","APOB","APOE","SORT1","HMGCR","NPC1L1","ABCG8",
              "LPA","CETP","LDLRAP1","MYLIP","ANGPTL3"),
    hdl   = c("CETP","LIPC","LIPG","ABCA1","APOA1","APOA5","LCAT","SCARB1",
              "PLTP","LPL","ANGPTL4","CUBN"),
    logtg = c("LPL","APOA5","APOC3","ANGPTL3","ANGPTL4","ANGPTL8","GCKR",
              "LMF1","GPIHBP1","APOE","TRIB1","MLXIPL"),
    tc    = c("LDLR","APOE","APOB","PCSK9","HMGCR","CETP","ABCG8","NPC1L1",
              "SORT1","TM6SF2","LPA","TRIB1"),
    character(0))
}
## every path above can still be overridden individually
GENES_OUT  <- Sys.getenv("MAGMA_GENES_OUT", GENES_OUT)
WY_CSV     <- Sys.getenv("WY_CSV",          WY_CSV)
SETS       <- Sys.getenv("PATHWAY_SETS",    SETS)
GENE_LOC   <- Sys.getenv("GENE_LOC",        GENE_LOC)
PW_NAMES_F <- Sys.getenv("PATHWAY_NAMES",   PW_NAMES_F)
GENES_ANNOT <- Sys.getenv("MAGMA_GENES_ANNOT", GENES_ANNOT)

## optional: SNP-level GWAS results.  If absent, S_GWAS is set to 0 for every
## gene and the score reduces to the MAGMA + pathway layers.  Expected columns
## are configurable below.
GWAS_FILE <- Sys.getenv("GWAS_FILE", if (SPECIES == "human") {
  cand <- c(file.path(DATA_D, sprintf("lipids/%s_eur.neglog10.txt", TRAIT)),
            file.path(DATA_D, "ldl/ldl_eur.neglog10.txt"),
            file.path(DATA_D, sprintf("lipids/%s_eur.magma.pval", TRAIT)),
            file.path(DATA_D, "ldl/ldl_eur.magma.pval"))
  cand <- if (TRAIT == "ldl") cand else cand[c(1, 3)]
  hit <- cand[file.exists(cand)]
  if (length(hit)) hit[1] else ""
} else "")
CHUNK <- as.integer(Sys.getenv("ANNOT_CHUNK", "2000"))
GWAS_SNP_COL <- Sys.getenv("GWAS_SNP_COL", "SNP")
GWAS_CHR_COL <- Sys.getenv("GWAS_CHR_COL", "Chr")
GWAS_POS_COL <- Sys.getenv("GWAS_POS_COL", "Pos")
GWAS_P_COL   <- Sys.getenv("GWAS_P_COL",   "p")
WINDOW       <- as.integer(Sys.getenv("GWAS_WINDOW", "10000"))

## pathway statistic to rank on (smaller = stronger): z_std (size-standardised
## WY statistic, the default), log10_p_omnibus (raw omnibus), or a p-value column
PCOL <- Sys.getenv("PATHWAY_P_COL", "z_std")

OUT_DIR <- file.path(BASE, "catfish2/results")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

## ---------------------------------------------------------------- helpers ---
## rank 1 = smallest value = strongest evidence; NA stays NA
rank_of <- function(p) {
  r <- rep(NA_real_, length(p))
  ok <- !is.na(p)
  r[ok] <- rank(p[ok], ties.method = "average")
  r
}
## -log10(rank/N), in [0, log10(N)]; genes with no evidence of this kind get 0
rank_score <- function(p) {
  out <- -log10(rank_of(p) / sum(!is.na(p)))
  out[is.na(out)] <- 0
  out
}
need <- function(f, what) if (!file.exists(f))
  stop(sprintf("%s not found: %s", what, f), call. = FALSE)

## ------------------------------------------------------------ gene layer ---
need(GENES_OUT, "MAGMA .genes.out")
g <- data.table::fread(GENES_OUT, data.table = FALSE)
pcol_magma <- if ("P_MULTI" %in% names(g)) "P_MULTI" else "P"
idcol <- if ("GENE" %in% names(g)) "GENE" else names(g)[1]
genes <- data.frame(gene = as.character(g[[idcol]]),
                    p_magma = as.numeric(g[[pcol_magma]]),
                    stringsAsFactors = FALSE)
genes <- genes[!is.na(genes$p_magma), ]
genes <- genes[!duplicated(genes$gene), ]
cat(sprintf("MAGMA genes: %d (p column %s)\n", nrow(genes), pcol_magma))
genes$S_MAGMA    <- rank_score(genes$p_magma)
genes$rank_magma <- rank_of(genes$p_magma)
N_MAGMA <- sum(!is.na(genes$p_magma))

## --------------------------------------------------------- pathway layer ---
need(WY_CSV, "Westfall-Young results"); need(SETS, "pathway set file")
wy <- data.table::fread(WY_CSV, data.table = FALSE)
if (!PCOL %in% names(wy))
  stop(sprintf("column %s not in %s (have: %s)", PCOL, WY_CSV,
               paste(names(wy), collapse = ", ")), call. = FALSE)
wy <- wy[!is.na(wy[[PCOL]]), ]
pw_p <- stats::setNames(as.numeric(wy[[PCOL]]), as.character(wy$pathway_id))

## sets file: one pathway per line, "<id> <gene> <gene> ..."
lines <- readLines(SETS, warn = FALSE)
lines <- lines[nzchar(trimws(lines))]
members <- lapply(strsplit(trimws(lines), "[ \t]+"), function(x) x)
pw_id   <- vapply(members, `[`, character(1), 1L)
pw_gene <- lapply(members, function(x) x[-1L])
names(pw_gene) <- pw_id
pw_gene <- pw_gene[names(pw_gene) %in% names(pw_p)]
cat(sprintf("pathways with a %s value and members: %d\n", PCOL, length(pw_gene)))

## best (smallest) pathway p per gene, number of pathways per gene
long <- data.frame(
  gene    = unlist(pw_gene, use.names = FALSE),
  pathway = rep(names(pw_gene), lengths(pw_gene)),
  stringsAsFactors = FALSE)
long$p <- pw_p[long$pathway]
agg <- stats::aggregate(p ~ gene, data = long, FUN = min)
npw <- stats::aggregate(pathway ~ gene, data = long,
                        FUN = function(x) length(unique(x)))
names(npw)[2] <- "n_pathways"
best <- merge(agg, npw, by = "gene", all = TRUE)
names(best)[names(best) == "p"] <- "p_pathway"
best$best_pathway <- long$pathway[match(
  paste(best$gene, best$p_pathway), paste(long$gene, long$p))]

genes <- merge(genes, best, by = "gene", all.x = TRUE)
genes$n_pathways[is.na(genes$n_pathways)] <- 0L
## rank-score pathway evidence only among annotated genes, then zero the rest
ann <- !is.na(genes$p_pathway)
genes$S_PATH <- 0
genes$S_PATH[ann] <- rank_score(genes$p_pathway[ann])
genes$rank_pathway <- NA_real_
genes$rank_pathway[ann] <- rank_of(genes$p_pathway[ann])
N_PATH <- sum(ann)
genes$B_hub <- log10(genes$n_pathways + 1)

## ----------------------------------------------------------- GWAS layer ---
genes$S_GWAS <- 0
genes$p_gwas <- NA_real_
genes$rank_gwas <- NA_real_
N_GWAS <- 0L

## Preferred route: MAGMA's .genes.annot already maps each gene to the SNPs it
## was given, using the same window as the gene analysis. That is both exactly
## consistent with the gene-level layer and the only option when the GWAS file
## carries no positions -- MAGMA's --pval format is SNP/P/N only.
## Fallback: a GWAS file with chromosome and position, mapped by window here.
if (nzchar(GWAS_FILE) && file.exists(GWAS_FILE) &&
    nzchar(GENES_ANNOT) && file.exists(GENES_ANNOT)) {
  cat("GWAS layer: using MAGMA gene annotation\n")
  pv <- data.table::fread(GWAS_FILE, data.table = TRUE, showProgress = FALSE)
  snp_col <- intersect(c("SNP","snp","ID","rsid"), names(pv))[1]
  ## Prefer -log10(p) when it is available. GLGC's strongest variants reach
  ## -log10(p) = 1725, far below the smallest representable double, so a
  ## probability column is floored and the strongest genes tie. Larger
  ## -log10(p) means stronger, so negate it to keep "smaller is stronger".
  nl_col <- intersect(c("NEGLOG10P","neglog10p","pvalue_neg_log10","NEG_LOG10_P"),
                      names(pv))[1]
  p_col  <- intersect(c("P","p","PVAL","pval","P_VALUE"), names(pv))[1]
  if (is.na(snp_col) || (is.na(p_col) && is.na(nl_col)))
    stop("could not find SNP and a p-value column in ", GWAS_FILE, call. = FALSE)
  if (!is.na(nl_col)) {
    cat("  ranking on -log10(p) (unfloored)\n")
    pv <- pv[, c(snp_col, nl_col), with = FALSE]
    data.table::setnames(pv, c("SNP", "P"))
    pv[, P := -P]
  } else {
    cat("  ranking on the probability column -- check for floor ties\n")
    pv <- pv[, c(snp_col, p_col), with = FALSE]
    data.table::setnames(pv, c("SNP", "P"))
  }
  pv <- pv[is.finite(P)]
  data.table::setkeyv(pv, "SNP")
  cat(sprintf("  %d SNPs with a p-value\n", nrow(pv)))

  ## stream the annotation in blocks so memory stays bounded: the long
  ## gene x SNP table for a whole genome would be tens of millions of rows
  con <- file(GENES_ANNOT, "r"); on.exit(close(con), add = TRUE)
  acc <- list(); nline <- 0L
  repeat {
    ln <- readLines(con, n = CHUNK, warn = FALSE)
    if (!length(ln)) break
    ln <- ln[!startsWith(ln, "#")]
    if (!length(ln)) next
    nline <- nline + length(ln)
    parts <- strsplit(ln, "[ \t]+")
    gids  <- vapply(parts, function(x) x[1], character(1))
    snps  <- lapply(parts, function(x) if (length(x) > 2L) x[-(1:2)] else character(0))
    keep  <- lengths(snps) > 0L
    if (!any(keep)) next
    dt <- data.table::data.table(
      gene = rep(gids[keep], lengths(snps[keep])),
      SNP  = unlist(snps[keep], use.names = FALSE))
    dt <- pv[dt, on = "SNP", nomatch = 0L]
    if (nrow(dt)) acc[[length(acc) + 1L]] <- dt[, list(p = min(P)), by = gene]
    cat(sprintf("\r  genes scanned: %d", nline))
  }
  cat("\n")
  if (length(acc)) {
    allg <- data.table::rbindlist(acc)
    allg <- allg[, list(p = min(p)), by = gene]
    genes$p_gwas <- allg$p[match(genes$gene, allg$gene)]
  }
} else if (nzchar(GWAS_FILE) && file.exists(GWAS_FILE) && file.exists(GENE_LOC)) {
  cat("GWAS layer: mapping SNPs to genes by position window\n")
  loc <- data.table::fread(GENE_LOC, header = FALSE, data.table = FALSE)
  names(loc)[1:4] <- c("gene", "chr", "start", "end")
  loc$chr <- as.character(loc$chr)
  gw <- data.table::fread(GWAS_FILE, data.table = FALSE)
  gw <- data.frame(chr = as.character(gw[[GWAS_CHR_COL]]),
                   pos = as.numeric(gw[[GWAS_POS_COL]]),
                   p   = as.numeric(gw[[GWAS_P_COL]]),
                   stringsAsFactors = FALSE)
  gw <- gw[!is.na(gw$p), ]
  best_snp <- rep(NA_real_, nrow(loc))
  for (cc in unique(loc$chr)) {
    li <- which(loc$chr == cc); gi <- gw[gw$chr == cc, ]
    if (!nrow(gi) || !length(li)) next
    gi <- gi[order(gi$pos), ]
    for (i in li) {
      sel <- gi$pos >= (loc$start[i] - WINDOW) & gi$pos <= (loc$end[i] + WINDOW)
      if (any(sel)) best_snp[i] <- min(gi$p[sel])
    }
  }
  loc$p_gwas <- best_snp
  genes$p_gwas <- loc$p_gwas[match(genes$gene, loc$gene)]
} else {
  cat("GWAS layer: no usable SNP file -- S_GWAS = 0 for all genes\n")
  cat("  set GWAS_FILE (and MAGMA_GENES_ANNOT, or a file with CHR/POS columns)\n")
}

hit <- !is.na(genes$p_gwas)
if (any(hit)) {
  genes$S_GWAS[hit]    <- rank_score(genes$p_gwas[hit])
  genes$rank_gwas[hit] <- rank_of(genes$p_gwas[hit])
  N_GWAS <- sum(hit)
  cat(sprintf("  %d genes with an overlapping SNP\n", N_GWAS))
  mn <- min(genes$p_gwas[hit])
  nt <- sum(genes$p_gwas[hit] <= mn * (1 + 1e-9))
  if (nt > 1L)
    cat(sprintf("  WARNING: %d genes tie at the strongest value (%.3g) -- the\n          GWAS layer cannot order them\n", nt, mn))
}

## ------------------------------------------------------------- composite ---
genes$Score <- genes$S_GWAS + genes$S_MAGMA + 0.5 * genes$S_PATH +
               0.3 * genes$B_hub
genes$layers <- (genes$S_GWAS > 0) + (genes$S_MAGMA > 0) + (genes$S_PATH > 0)
genes <- genes[order(-genes$Score), ]
genes$rank <- seq_len(nrow(genes))

## gene symbols (human Entrez IDs are unreadable on their own)
if (!is.na(SYMBOL_COL) && file.exists(GENE_LOC)) {
  gl <- data.table::fread(GENE_LOC, header = FALSE, data.table = FALSE)
  if (ncol(gl) >= SYMBOL_COL) {
    sym <- stats::setNames(as.character(gl[[SYMBOL_COL]]), as.character(gl[[1]]))
    genes$symbol <- unname(sym[genes$gene])
  }
}

## readable pathway names
if (nzchar(PW_NAMES_F) && file.exists(PW_NAMES_F)) {
  nm <- utils::read.delim(PW_NAMES_F, stringsAsFactors = FALSE,
                          check.names = FALSE)
  if (all(c("Pathway-id", "Pathway-name") %in% names(nm))) {
    nm <- unique(nm[, c("Pathway-id", "Pathway-name")])
    key <- stats::setNames(nm[["Pathway-name"]], nm[["Pathway-id"]])
    genes$best_pathway_name <- unname(key[genes$best_pathway])
  }
}

out <- file.path(OUT_DIR, sprintf("06_candidate_genes_%s.csv", TAG))
utils::write.csv(genes, out, row.names = FALSE)
cat(sprintf("\nwrote %s  (%d genes)\n", out, nrow(genes)))

cat(sprintf("\nlayer sizes: MAGMA %d genes | pathway-annotated %d | GWAS-overlapping %d\n",
            N_MAGMA, N_PATH, N_GWAS))
cols <- intersect(c("rank","gene","symbol","Score","rank_gwas","rank_magma",
                    "rank_pathway","p_magma","p_pathway","best_pathway",
                    "best_pathway_name","n_pathways","layers"),
                  names(genes))
cat("\nTop 20 candidate genes:\n")
print(utils::head(genes[, cols], 20), row.names = FALSE)

## positive controls: genes the literature already implicates
lookup <- if ("symbol" %in% names(genes)) genes$symbol else genes$gene
TARGETS <- unlist(strsplit(Sys.getenv("TARGET_GENES",
                                      paste(TARGETS, collapse = ",")), ","))
cat("\nPositive-control genes:\n")
hits <- 0L
for (t in TARGETS) {
  i <- which(lookup == t)
  if (!length(i)) { cat(sprintf("  %-12s not in the MAGMA gene table\n", t)); next }
  i <- i[1]; hits <- hits + 1L
  fmt <- function(v, n) if (is.na(v)) "      -" else sprintf("%6.0f/%d", v, n)
  cat(sprintf("  %-12s overall %6d of %d | GWAS %s | MAGMA %s | pathway %s | %s\n",
              t, genes$rank[i], nrow(genes),
              fmt(genes$rank_gwas[i], N_GWAS),
              fmt(genes$rank_magma[i], N_MAGMA),
              fmt(genes$rank_pathway[i], N_PATH),
              ifelse(is.na(genes$best_pathway[i]), "no pathway",
                     substr(genes$best_pathway[i], 1, 40))))
}
if (hits > 1L) {
  idx <- match(TARGETS, lookup); idx <- idx[!is.na(idx)]
  rk <- genes$rank[idx]
  cat(sprintf("\n  %d controls found: median rank %.0f of %d (random expectation %.0f)\n",
              length(rk), stats::median(rk), nrow(genes), nrow(genes) / 2))
  cat(sprintf("  in the top 1%%: %d of %d\n",
              sum(rk <= nrow(genes) * 0.01), length(rk)))
}
