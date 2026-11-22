#!/usr/bin/env Rscript
## =============================================================================
## 06_logo_bap.R <trait> -- leave-one-gene-out pathway evidence for the BAP stem
## composition gene score (Section 3.10), as scripts/hpc/gc_real/wy_logo.R does
## for the lipid controls. For each gene g to check and each tested pathway P
## that contains it, P and P without g are scored together in one
## Westfall-Young run, so the change in z_std from removing g is measured on
## the same null draws. 47_bap_logo.R then gives g the pathway evidence of its
## pathways without itself and recomputes its rank in the gene score.
##
## Genes checked: every gene with pathway evidence among the top LOGO_TOPN
## (default 100) of the trait's gene score (prioritize/<trait>_gene_prioritize.csv),
## which includes SbALDH3E1 (SORBI_3004G300800) and the thioredoxin reductase
## (SORBI_3004G300900).
## Inputs and settings are those of 03_catfish_bap.R: MAGMA SNP-wise mean
## .genes.out (10 chromosomes), gene-level GC, catfish/<trait>/gene_cor_pairs.txt,
## SorghumBicolorCyc (.SORBI) pathways, B = 10000, 5 components, tau {0.01,
## 0.05, 0.5}, seed 20260706; pathways: those in catfish/<trait>/<trait>_wy.csv;
## the reduced sets may have one gene fewer than the 5-gene minimum.
## DEV_WY as in 03 (default: the gc_real dev copy, if present).
##   Rscript 06_logo_bap.R Brix_maturity
## Writes logo/<trait>_logo.csv: pathway_id, base_pathway, removed_gene, z_std, n_genes
## =============================================================================
suppressPackageStartupMessages(library(CATFISH))
TRAIT <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(TRAIT)) stop("usage: Rscript 06_logo_bap.R <trait>", call. = FALSE)
ROOT <- Sys.getenv("BAP_ROOT", "/share/maize/ntanduk/BAP_gwas")
PWF  <- Sys.getenv("PATHWAY_FILE",
          "/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/inputs/sorghumbicolorcyc_pathways.20230103.SORBI")
DEV  <- Sys.getenv("DEV_WY", "/share/maize/ntanduk/CATFISH/gc_real/dev/westfall_young.R")
B    <- as.integer(Sys.getenv("WY_B", "10000"))
TOPN <- as.integer(Sys.getenv("LOGO_TOPN", "100"))
SEED <- 20260706L; TAU <- c(0.01, 0.05, 0.5); MIN <- 5L; MAX <- 200L
CHROMS <- sprintf("%02d", 1:10)
if (nzchar(DEV) && file.exists(DEV)) {
  dev_env <- new.env(parent = asNamespace("CATFISH"))
  sys.source(DEV, envir = dev_env)
  catfish_wy_pathways <- dev_env$catfish_wy_pathways
  cat("using catfish_wy_pathways() from", DEV, "\n")
}

MD   <- file.path(ROOT, "magma", TRAIT)
CD   <- file.path(ROOT, "catfish", TRAIT)
OUT  <- file.path(ROOT, "logo"); dir.create(OUT, showWarnings = FALSE)
fo   <- file.path(MD, sprintf("%s_chr%s.snp_wise_mean.genes.out", TRAIT, CHROMS))
CORP <- file.path(CD, "gene_cor_pairs.txt")
WYF  <- file.path(CD, sprintf("%s_wy.csv", TRAIT))
PRF  <- file.path(ROOT, "prioritize", sprintf("%s_gene_prioritize.csv", TRAIT))
miss <- c(fo, CORP, WYF, PRF, PWF)[!file.exists(c(fo, CORP, WYF, PRF, PWF))]
if (length(miss)) stop("missing:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

## ---- gene statistics and GC, as 03 ---------------------------------------------------
g <- do.call(rbind, lapply(fo, utils::read.table, header = TRUE, stringsAsFactors = FALSE, comment.char = "#"))
g <- g[!duplicated(g$GENE) & is.finite(g$P), ]
z <- stats::qnorm(log(pmax(g$P, 1e-320)), lower.tail = FALSE, log.p = TRUE)
g$P_GC <- stats::pnorm(z - stats::median(z), lower.tail = FALSE)

## ---- pathways, as 03; tested = those in the paper's result -----------------------------
raw <- utils::read.delim(PWF, stringsAsFactors = FALSE, check.names = FALSE)
pw0 <- unique(data.frame(pathway_id = raw[["Pathway-id"]], gene_id = raw[["Gene-name"]]))
pw0 <- pw0[nzchar(pw0$gene_id) & pw0$gene_id != "unknown", ]
pw  <- split(pw0$gene_id, pw0$pathway_id)
tested <- utils::read.csv(WYF, stringsAsFactors = FALSE)$pathway_id

## ---- genes to check: pathway-supported genes in the top TOPN of the gene score ----------
pr <- utils::read.csv(PRF, stringsAsFactors = FALSE)
chk <- pr$gene[pr$rank <= TOPN & !is.na(pr$best_pathway) & nzchar(pr$best_pathway)]
chk <- unique(c(chk, intersect(c("SORBI_3004G300800", "SORBI_3004G300900"), pr$gene)))
sets <- list(); map <- list()
for (id in chk) for (P in tested[vapply(tested, function(p) id %in% pw[[p]], logical(1))]) {
  sets[[P]] <- pw[[P]]
  k <- paste0(P, "__minus__", id)
  sets[[k]] <- setdiff(pw[[P]], id)
  map[[k]] <- data.frame(pathway_id = c(P, k), base_pathway = P, removed_gene = c(NA, id),
                         stringsAsFactors = FALSE)
}
map <- unique(do.call(rbind, map))
cat(sprintf("[%s] %d genes checked in %d tested pathways -> %d sets, B = %d\n",
            TRAIT, length(chk), length(unique(map$base_pathway)), length(sets), B))

t0 <- Sys.time()
res <- catfish_wy_pathways(g, pathways = sets, gene_col = "GENE", p_col = "P_GC", chr_col = "CHR",
                           cor_file = CORP, B = B, use_tfisher = TRUE, tau_grid = TAU,
                           min_genes = MIN - 1L, max_genes = MAX, seed = SEED, verbose = FALSE)
out <- merge(map, res[, c("pathway_id", "z_std", "n_genes")], by = "pathway_id", all.x = TRUE)
fw <- file.path(OUT, sprintf("%s_logo.csv", TRAIT))
utils::write.csv(out, fw, row.names = FALSE)
cat(sprintf("[%s] done in %.1f min -> %s\n", TRAIT, as.numeric(difftime(Sys.time(), t0, units = "mins")), fw))
