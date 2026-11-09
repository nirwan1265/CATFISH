#!/usr/bin/env Rscript
## =============================================================================
## 25_sorghum_gene_prioritize.R -- candidate genes for the sorghum BAP traits
## with CATFISH::gene_prioritize(), and how the pre-specified established genes
## rank under the combined score versus the GWAS layer alone.
##
## Inputs (all from the runs already done; nothing is recomputed):
##   MAGMA     <MAGMA_ROOT>/<TRAIT>/<TRAIT>_chrNN.multi_snp_wise.genes.out
##             gene p column P_SNPWISE_MEAN (the column 23 used; see 24)
##   pathways  ../results/sorghum_bap/<TRAIT>_P_SNPWISE_MEAN_5comp_westfall_young.csv
##             (from 23) and the same SorghumBicolorCyc .SORBI set file
##   GWAS      <MAGMA_ROOT>/<TRAIT>/<TRAIT>_chrNN.pval.txt -- the SNP p-values
##             MAGMA itself read (GEMMA p_wald); SNPs with N < 50 dropped, as
##             MAGMA did. CHR and BP are parsed from the SNP ID (Chr01pos236.1).
##   SNP -> gene: by position, using the gene START/STOP that MAGMA reports in
##             .genes.out. Those coordinates already include the +/-25 kb window
##             of the MAGMA annotation (e.g. SORBI_3008G000100 is 10,044-14,451
##             in sorghum.genes.loc and 0-39,451 in .genes.out), so window = 0
##             here gives each gene exactly the SNPs MAGMA used (checked on
##             chr 8: SNP counts match NSNPS for 100% of genes).
##
## Usage (from catfish2/scripts):
##   Rscript 25_sorghum_gene_prioritize.R                     # both traits
##   TRAITS=Lignin_DM_perc Rscript 25_sorghum_gene_prioritize.R
## Env: MAGMA_ROOT, CATFISH_BASE, PATHWAY_FILE, CHROMS (quick checks only)
## Writes ../results/sorghum_bap/25_<TRAIT>_gene_prioritize.csv   (all genes)
##        ../results/sorghum_bap/25_established_genes.csv         (both traits)
## =============================================================================
suppressPackageStartupMessages({ library(CATFISH); library(data.table) })

MAGMA_ROOT <- Sys.getenv("MAGMA_ROOT", path.expand("~/Documents/Research/results/CATFISH/MAGMA"))
BASE       <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/CATFISH"))
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/sorghumbicolorcyc_pathways.20230103.SORBI"))
RES_DIR    <- Sys.getenv("RES_DIR", "../results/sorghum_bap")
TRAITS     <- strsplit(Sys.getenv("TRAITS", "Dry_tons_per_acre,Lignin_DM_perc"), ",")[[1]]
CHROMS     <- as.integer(strsplit(Sys.getenv("CHROMS", "1,2,3,4,5,6,7,8,9,10"), ",")[[1]])
P_COL      <- "P_SNPWISE_MEAN"
MIN_N      <- 50L                      # MAGMA discarded SNPs with N < 50

## established genes, named before any result was seen (same list as 24)
ESTABLISHED <- data.frame(
  trait = c(rep("Dry_tons_per_acre", 8), rep("Lignin_DM_perc", 3)),
  gene  = c("Ma1/PRR37", "Ma3/PHYB", "Ma5/PHYC", "Ma6/GHD7", "Ma2",
            "Dw1", "Dw2", "Dw3", "bmr2/4CL", "bmr6/CAD", "bmr12/COMT"),
  id    = c("SORBI_3006G057866", "SORBI_3001G394400", "SORBI_3001G087100",
            "SORBI_3006G004400", "SORBI_3002G302700", "SORBI_3009G229800",
            "SORBI_3006G067700", "SORBI_3007G163800", "SORBI_3004G062500",
            "SORBI_3004G071000", "SORBI_3007G047300"),
  stringsAsFactors = FALSE)

## pathway membership: exactly the sets 23 tested
pw_raw <- fread(PATHWAY_FILE, sep = "\t", quote = "", showProgress = FALSE)
pw <- unique(data.frame(pathway_id = pw_raw[["Pathway-id"]],
                        gene_id    = pw_raw[["Gene-name"]], stringsAsFactors = FALSE))
pw <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]

est_all <- list()
for (TRAIT in TRAITS) {
  t0 <- Sys.time()
  cat(sprintf("\n================ %s ================\n", TRAIT))
  tdir <- file.path(MAGMA_ROOT, TRAIT)
  gfiles <- file.path(tdir, sprintf("%s_chr%02d.multi_snp_wise.genes.out", TRAIT, CHROMS))
  pfiles <- file.path(tdir, sprintf("%s_chr%02d.pval.txt", TRAIT, CHROMS))
  wy_file <- file.path(RES_DIR, sprintf("%s_%s_5comp_westfall_young.csv", TRAIT, P_COL))
  miss <- c(gfiles, pfiles, wy_file)[!file.exists(c(gfiles, pfiles, wy_file))]
  if (length(miss)) stop("missing input files:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)

  mg <- rbindlist(lapply(gfiles, fread, showProgress = FALSE))
  mg <- mg[!duplicated(GENE)]
  loc <- data.frame(gene = mg$GENE, chr = as.character(mg$CHR),
                    start = mg$START, end = mg$STOP, stringsAsFactors = FALSE)

  gw <- rbindlist(lapply(pfiles, fread, showProgress = FALSE))
  gw <- gw[N >= MIN_N]
  gw[, `:=`(CHR = as.character(as.integer(sub("^Chr(\\d+)pos.*$", "\\1", SNP))),
            BP  = as.numeric(sub("^Chr\\d+pos(\\d+).*$", "\\1", SNP)))]
  cat(sprintf("GWAS: %s SNPs (N >= %d), min p %.3g\n",
              format(nrow(gw), big.mark = ","), MIN_N, min(gw$P)))

  cand <- gene_prioritize(
    magma = as.data.frame(mg), magma_gene_col = "GENE", magma_p_col = P_COL,
    pathway_results = wy_file, pathway_stat_col = "z_std",
    pathways = pw,
    gwas = as.data.frame(gw), gwas_snp_col = "SNP", gwas_p_col = "P",
    gwas_chr_col = "CHR", gwas_pos_col = "BP",
    gene_loc = loc, window = 0L,
    targets = ESTABLISHED$id[ESTABLISHED$trait == TRAIT],
    verbose = TRUE)

  ## readable pathway names
  wy <- fread(wy_file, showProgress = FALSE)
  cand$best_pathway_name <- wy$pathway_name[match(cand$best_pathway, wy$pathway_id)]
  out <- file.path(RES_DIR, sprintf("25_%s_gene_prioritize.csv", TRAIT))
  fwrite(cand, out)

  cat("\ntop 20 candidates:\n")
  print(head(data.frame(rank = cand$rank, gene = cand$gene, score = round(cand$score, 2),
                        gwas_rank = cand$gwas_rank, magma_rank = cand$magma_rank,
                        path_rank = cand$pathway_rank,
                        best_pathway = substr(cand$best_pathway_name, 1, 40)), 20),
        row.names = FALSE)

  ## established genes: combined vs single layers
  N <- nrow(cand)
  e <- ESTABLISHED[ESTABLISHED$trait == TRAIT, ]
  i <- match(e$id, cand$gene)
  e <- data.frame(e, n_genes = N,
                  combined_rank = cand$rank[i], gwas_rank = cand$gwas_rank[i],
                  magma_rank = cand$magma_rank[i], pathway_rank = cand$pathway_rank[i],
                  in_pathway = !is.na(cand$best_pathway[i]),
                  best_pathway = cand$best_pathway_name[i],
                  stringsAsFactors = FALSE)
  e$combined_top_pct <- round(100 * e$combined_rank / N, 1)
  e$gwas_top_pct     <- round(100 * e$gwas_rank / N, 1)
  e$magma_top_pct    <- round(100 * e$magma_rank / N, 1)
  est_all[[TRAIT]] <- e
  cat("\nestablished genes (rank out of", N, "genes; top x%):\n")
  print(e[, c("gene", "combined_rank", "gwas_rank", "magma_rank", "combined_top_pct",
              "gwas_top_pct", "magma_top_pct", "in_pathway")], row.names = FALSE)
  cat(sprintf("median rank: combined %.0f, GWAS alone %.0f, MAGMA alone %.0f\n",
              median(e$combined_rank, na.rm = TRUE), median(e$gwas_rank, na.rm = TRUE),
              median(e$magma_rank, na.rm = TRUE)))
  cat(sprintf("better under combined than GWAS alone: %d of %d\n",
              sum(e$combined_rank < e$gwas_rank, na.rm = TRUE), sum(!is.na(e$combined_rank))))
  cat(sprintf("wrote %s (%.1f min)\n", out, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
est <- do.call(rbind, est_all)
fwrite(est, file.path(RES_DIR, "25_established_genes.csv"))
cat(sprintf("\nwrote %s\n", file.path(RES_DIR, "25_established_genes.csv")))
