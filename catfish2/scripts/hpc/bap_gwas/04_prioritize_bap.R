#!/usr/bin/env Rscript
## =============================================================================
## 04_prioritize_bap.R -- CATFISH gene prioritisation for BAP traits, and the
## genes shared between them.
## Per trait, as catfish2/scripts/25_sorghum_gene_prioritize.R (the paper's
## sorghum gene score):
##   MAGMA     magma/<trait>/<trait>_chrNN.snp_wise_mean.genes.out (column P)
##   pathways  catfish/<trait>/<trait>_wy.csv, pathway statistic z_std, with the
##             SorghumBicolorCyc membership (.SORBI)
##   GWAS      the trait's GEMMA .assoc.txt (rs, chr, ps, p_wald); SNPs mapped to
##             genes by position with the START/STOP that MAGMA reports, which
##             already include the +/- 25 kb window (window = 0 here)
##   score     CATFISH::gene_prioritize(), weights GWAS 1, MAGMA 1, pathway W_PATH
##             (default 1), hub W_HUB (default 0); until 2026-11-25 0.5 and 0.3
## Shared genes: a gene's rank in every trait; genes in the top K of at least
## MIN_TRAITS traits (K = 50, 100, 200).
## Writes prioritize/<trait>_gene_prioritize.csv, prioritize/shared_genes.csv
##   Rscript 04_prioritize_bap.R
## Env: TRAITS (comma-separated; default the five below), BAP_ROOT, RUNS_DIR,
##      MAGMA_DIR, CATFISH_DIR (default runs/, magma/, catfish/ under BAP_ROOT),
##      PATHWAY_FILE, MIN_TRAITS (default 3), SHARED_NAME (file name of the
##      shared-gene table; default shared_genes.csv)
## =============================================================================
suppressPackageStartupMessages({ library(CATFISH); library(data.table) })
ROOT  <- Sys.getenv("BAP_ROOT", "/share/maize/ntanduk/BAP_gwas")
RUNS  <- Sys.getenv("RUNS_DIR", file.path(ROOT, "runs"))
MAGD  <- Sys.getenv("MAGMA_DIR", file.path(ROOT, "magma"))
CATD  <- Sys.getenv("CATFISH_DIR", file.path(ROOT, "catfish"))
OUT   <- file.path(ROOT, "prioritize"); dir.create(OUT, showWarnings = FALSE)
PWF   <- Sys.getenv("PATHWAY_FILE",
          "/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy/inputs/sorghumbicolorcyc_pathways.20230103.SORBI")
TRAITS <- strsplit(Sys.getenv("TRAITS", "ADF_pct_DM,NDF_pct_DM,NFC_pct_DM,Brix_milk,Brix_maturity"), ",")[[1]]
MIN_TRAITS <- as.integer(Sys.getenv("MIN_TRAITS", "3"))
SHARED_NAME <- Sys.getenv("SHARED_NAME", "shared_genes.csv")
KS <- c(50, 100, 200)
CHROMS <- sprintf("%02d", 1:10)

raw <- fread(PWF, sep = "\t", quote = "", showProgress = FALSE)
pw  <- unique(data.frame(pathway_id = raw[["Pathway-id"]], gene_id = raw[["Gene-name"]], stringsAsFactors = FALSE))
pw  <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]
pw_name <- tapply(raw[["Pathway-name"]], raw[["Pathway-id"]], function(x) x[1])

ranks <- list()
for (TRAIT in TRAITS) {
  t0 <- Sys.time(); cat(sprintf("\n==== %s ====\n", TRAIT))
  fo  <- file.path(MAGD, TRAIT, sprintf("%s_chr%s.snp_wise_mean.genes.out", TRAIT, CHROMS))
  wyf <- file.path(CATD, TRAIT, sprintf("%s_wy.csv", TRAIT))
  asc <- list.files(file.path(RUNS, TRAIT), pattern = "\\.assoc\\.txt$", recursive = TRUE, full.names = TRUE)
  miss <- c(fo, wyf)[!file.exists(c(fo, wyf))]
  if (length(miss) || length(asc) != 1) {
    cat("skip: missing", paste(miss, collapse = ", "), if (length(asc) != 1) sprintf("(%d .assoc.txt)", length(asc)), "\n"); next }

  mg <- rbindlist(lapply(fo, fread, showProgress = FALSE)); mg <- mg[!duplicated(GENE) & is.finite(P)]
  loc <- data.frame(gene = mg$GENE, chr = as.character(as.integer(mg$CHR)),
                    start = mg$START, end = mg$STOP, stringsAsFactors = FALSE)
  gw <- fread(asc, select = c("chr", "rs", "ps", "p_wald"), showProgress = FALSE)
  gw <- data.frame(SNP = gw$rs, CHR = as.character(as.integer(gw$chr)), BP = gw$ps, P = gw$p_wald)
  cat(sprintf("MAGMA %d genes; GWAS %s SNPs, min p %.3g\n", nrow(mg), format(nrow(gw), big.mark = ","), min(gw$P)))

  cand <- gene_prioritize(
    magma = as.data.frame(mg), magma_gene_col = "GENE", magma_p_col = "P",
    pathway_results = wyf, pathway_stat_col = "z_std", pathways = pw,
    gwas = gw, gwas_snp_col = "SNP", gwas_p_col = "P", gwas_chr_col = "CHR", gwas_pos_col = "BP",
    gene_loc = loc, window = 0L,
    w_pathway = as.numeric(Sys.getenv("W_PATH", "1")), w_hub = as.numeric(Sys.getenv("W_HUB", "0")), verbose = FALSE)
  cand$best_pathway_name <- unname(pw_name[cand$best_pathway])
  cand$chr <- loc$chr[match(cand$gene, loc$gene)]
  cand$start <- loc$start[match(cand$gene, loc$gene)]
  fwrite(cand, file.path(OUT, sprintf("%s_gene_prioritize.csv", TRAIT)))
  ranks[[TRAIT]] <- data.table(gene = cand$gene, trait = TRAIT, rank = cand$rank,
                               gwas_rank = cand$gwas_rank, magma_rank = cand$magma_rank,
                               pathway_rank = cand$pathway_rank, best_pathway = cand$best_pathway_name,
                               chr = cand$chr, start = cand$start)
  cat("top 10:\n")
  print(head(data.frame(rank = cand$rank, gene = cand$gene, gwas = cand$gwas_rank, magma = cand$magma_rank,
                        path = cand$pathway_rank, best_pathway = substr(cand$best_pathway_name, 1, 40)), 10),
        row.names = FALSE)
  cat(sprintf("%.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

## ---- genes shared across traits ------------------------------------------------
R <- rbindlist(ranks)
if (!nrow(R)) stop("no trait finished", call. = FALSE)
W <- dcast(R, gene + chr + start ~ trait, value.var = "rank")
nt <- length(unique(R$trait))
for (k in KS) W[[sprintf("n_top%d", k)]] <- rowSums(W[, unique(R$trait), with = FALSE] <= k, na.rm = TRUE)
bp <- R[order(rank), .(best_pathway = best_pathway[!is.na(best_pathway)][1]), by = gene]
W <- merge(W, bp, by = "gene", all.x = TRUE)
W[, mean_rank := rowMeans(.SD, na.rm = TRUE), .SDcols = unique(R$trait)]
setorder(W, -n_top100, mean_rank)
fwrite(W, file.path(OUT, SHARED_NAME))
options(width = 200)
cat(sprintf("\n==== genes in the top 100 of at least %d of %d traits ====\n", MIN_TRAITS, nt))
print(W[n_top100 >= MIN_TRAITS], row.names = FALSE)
cat(sprintf("\ngenes in the top 50 / 100 / 200 of all %d traits: %d / %d / %d\n", nt,
            sum(W$n_top50 == nt), sum(W$n_top100 == nt), sum(W$n_top200 == nt)))
cat("wrote", file.path(OUT, SHARED_NAME), "\n")
