#!/usr/bin/env Rscript
## =============================================================================
## 31_fly_gene_prioritize.R -- candidate genes for the DGRP body-weight and
## triglyceride traits with CATFISH::gene_prioritize(), and how pre-specified
## established genes rank under the combined score versus each layer alone.
## Same analysis as 25_sorghum_gene_prioritize.R, on the fly inputs:
##   MAGMA     <MAGMA_DIR>/<TRAIT>.genes.out, column P (SNP-wise mean)
##   pathways  ../results/dgrp/<TRAIT>_reactome_5comp_westfall_young.csv (from 30)
##             and the same Reactome fly set file
##   GWAS      <GWAS_DIR>/<TRAIT>.assoc.txt (vcf2gwas/GEMMA p_wald); N per SNP =
##             lines with the phenotype minus n_miss, and SNPs with N < 50 are
##             dropped, as MAGMA did ("discarding SNPs with N < 50" in its log)
##   SNP -> gene: by position, using the gene START/STOP in .genes.out, which
##             already include the +/-1 kb MAGMA window, so window = 0 gives each
##             gene the SNPs MAGMA used
## Usage (from catfish2/scripts):
##   Rscript 31_fly_gene_prioritize.R
##   TRAITS=TG_control Rscript 31_fly_gene_prioritize.R
##   TRAITS=BW_serrano_female_web,BW_serrano_male_web,TG_serrano_web Rscript 31_fly_gene_prioritize.R
##   (web-server traits: N from gwas/output/<TRAIT>.nobs; TG_* traits use the TG
##    established genes, all others the body-weight genes)
## Writes ../results/dgrp/31_<TRAIT>_gene_prioritize.csv  (all genes)
##        ../results/dgrp/31_established_genes.csv         (all traits)
## =============================================================================
suppressPackageStartupMessages({ library(CATFISH); library(data.table) })
BASE      <- Sys.getenv("CATFISH_BASE", path.expand("~/Documents/Github/MAGCAT"))
MAGMA_DIR <- Sys.getenv("MAGMA_DIR", file.path(BASE, "dgrp_bw_tg/magma"))
GWAS_DIR  <- Sys.getenv("GWAS_DIR",  file.path(BASE, "dgrp_bw_tg/gwas/output"))
PHENO_TSV <- Sys.getenv("PHENO_TSV", file.path(BASE, "dgrp_bw_tg/geno/pheno.tsv"))
PATHWAY_FILE <- Sys.getenv("PATHWAY_FILE",
                  file.path(BASE, "inst/extdata/pathway/dme_reactome_all_levels.tsv"))
GENE_MAP  <- file.path(BASE, "inst/extdata/dmel.flybase.gene_map.tsv")
RES_DIR   <- "../results/dgrp"
TRAITS    <- strsplit(Sys.getenv("TRAITS", "BW_control_female,BW_control_male,TG_control"), ",")[[1]]
MIN_N     <- 50L

## established genes, named before any gene-level result was seen.
## Body weight: insulin/TOR growth pathway (size control in Drosophila).
## Triglycerides: lipid storage and lipolysis, plus the insulin pathway core.
BW_GENES <- c("InR", "chico", "Pi3K92E", "Akt", "foxo", "S6k", "Thor", "Ilp2",
              "Pten", "gig", "Lnk", "Myc")
TG_GENES <- c("bmm", "Lsd-1", "Lsd-2", "AkhR", "Akh", "mdy", "Lpin", "adp",
              "SREBP", "InR", "chico", "foxo")
gm <- fread(GENE_MAP, showProgress = FALSE)
sym2id <- setNames(gm$FBgn, gm$symbol); id2sym <- setNames(gm$symbol, gm$FBgn)
ESTABLISHED <- rbindlist(lapply(TRAITS, function(t) {
  s <- if (grepl("^TG_", t)) TG_GENES else BW_GENES
  data.table(trait = t, gene = s, id = unname(sym2id[s]))
}))
if (anyNA(ESTABLISHED$id)) stop("no FBgn for: ", paste(unique(ESTABLISHED$gene[is.na(ESTABLISHED$id)]), collapse = ", "))

pw_raw <- fread(PATHWAY_FILE, sep = "\t", header = FALSE, quote = "", showProgress = FALSE,
                col.names = c("gene_id", "pathway_id", "url", "name", "evidence", "species"))
pw <- unique(as.data.frame(pw_raw[, .(pathway_id, gene_id)]))
ph <- fread(PHENO_TSV, showProgress = FALSE)

est_all <- list()
for (TRAIT in TRAITS) {
  t0 <- Sys.time()
  cat(sprintf("\n================ %s ================\n", TRAIT))
  gfile  <- file.path(MAGMA_DIR, paste0(TRAIT, ".genes.out"))
  afile  <- file.path(GWAS_DIR, paste0(TRAIT, ".assoc.txt"))
  wy_file <- file.path(RES_DIR, sprintf("%s_reactome_5comp_westfall_young.csv", TRAIT))
  miss <- c(gfile, afile, wy_file)[!file.exists(c(gfile, afile, wy_file))]
  if (length(miss)) stop("missing input files:\n  ", paste(miss, collapse = "\n  "), call. = FALSE)
  ## lines in the GWAS: the trait's column of pheno.tsv (local vcf2gwas runs), or
  ## <TRAIT>.nobs written by 00_web_to_assoc.py (DGRP web-server runs, *_web)
  nobs_file <- file.path(GWAS_DIR, paste0(TRAIT, ".nobs"))
  if (TRAIT %in% names(ph)) {
    NOBS <- sum(!is.na(ph[[TRAIT]]))
  } else if (file.exists(nobs_file)) {
    NOBS <- as.integer(readLines(nobs_file, n = 1))
  } else stop(TRAIT, " is not a column of ", PHENO_TSV, " and has no ", nobs_file)

  mg <- fread(gfile, showProgress = FALSE)[!duplicated(GENE)]
  loc <- data.frame(gene = mg$GENE, chr = as.character(mg$CHR),
                    start = mg$START, end = mg$STOP, stringsAsFactors = FALSE)
  gw <- fread(afile, select = c("chr", "rs", "ps", "n_miss", "p_wald"), showProgress = FALSE)
  gw <- gw[!is.na(p_wald) & (NOBS - n_miss) >= MIN_N]
  gw <- data.frame(SNP = gw$rs, P = gw$p_wald, CHR = as.character(gw$chr), BP = gw$ps)
  cat(sprintf("GWAS: %s SNPs (N >= %d of %d lines), min p %.3g\n",
              format(nrow(gw), big.mark = ","), MIN_N, NOBS, min(gw$P)))

  cand <- gene_prioritize(
    magma = as.data.frame(mg), magma_gene_col = "GENE", magma_p_col = "P",
    pathway_results = wy_file, pathway_stat_col = "z_std",
    pathways = pw,
    gwas = gw, gwas_snp_col = "SNP", gwas_p_col = "P",
    gwas_chr_col = "CHR", gwas_pos_col = "BP",
    gene_loc = loc, window = 0L,
    targets = ESTABLISHED$id[ESTABLISHED$trait == TRAIT],
    verbose = TRUE)
  wy <- fread(wy_file, showProgress = FALSE)
  cand$symbol <- unname(id2sym[cand$gene])
  cand$best_pathway_name <- wy$pathway_name[match(cand$best_pathway, wy$pathway_id)]
  out <- file.path(RES_DIR, sprintf("31_%s_gene_prioritize.csv", TRAIT))
  fwrite(cand, out)
  cat("\ntop 20 candidates:\n")
  print(head(data.frame(rank = cand$rank, gene = cand$symbol, id = cand$gene,
                        score = round(cand$score, 2), gwas_rank = cand$gwas_rank,
                        magma_rank = cand$magma_rank, path_rank = cand$pathway_rank,
                        best_pathway = substr(cand$best_pathway_name, 1, 40)), 20),
        row.names = FALSE)

  N <- nrow(cand)
  e <- as.data.frame(ESTABLISHED[trait == TRAIT])
  i <- match(e$id, cand$gene)
  e <- data.frame(e, n_genes = N,
                  combined_rank = cand$rank[i], gwas_rank = cand$gwas_rank[i],
                  magma_rank = cand$magma_rank[i], pathway_rank = cand$pathway_rank[i],
                  in_pathway = !is.na(cand$best_pathway[i]),
                  best_pathway = cand$best_pathway_name[i], stringsAsFactors = FALSE)
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
## one file per run, named after its traits, so a serrano run does not
## overwrite the control table
est_file <- file.path(RES_DIR, if (identical(TRAITS, c("BW_control_female", "BW_control_male", "TG_control")))
  "31_established_genes.csv" else sprintf("31_established_genes_%s.csv", paste(TRAITS, collapse = "+")))
fwrite(est, est_file)
cat(sprintf("\nwrote %s\n", est_file))
