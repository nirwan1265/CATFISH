#!/usr/bin/env Rscript
## =============================================================================
## 32_artp_prep.R -- sARTP inputs for one lipid set (e.g. ldl_eur):
##   pathways   the Reactome pathways CATFISH tested for this set: genes of the
##              set's MAGMA .genes.out with finite P and ZSTAT, 10-200 of them
##              per pathway (as wy_gc_real.R selects them)
##   SNP->gene  MAGMA's g1000_<panel>.genes.annot (the annotation the lipid MAGMA
##              run used), restricted to the genes of those pathways
##   summary    the GLGC file MAGMA's p-values came from: rsID, REF, ALT,
##              EFFECT_SIZE (per ALT allele), SE, and p from pvalue_neg_log10
##              (the column MAGMA was given; the plain pvalue column holds values
##              below 1e-308 that do not parse as numbers), floored at 1e-300 as
##              for MAGMA, for SNPs in the panel and in the SNP->gene map
## Writes ARTP_IN/<ta>.summary.rds, <ta>.pathways.rds, <ta>.settings.tsv.
##   Rscript 32_artp_prep.R ldl_eur
## =============================================================================
suppressPackageStartupMessages(library(data.table))
TA <- commandArgs(trailingOnly = TRUE)[1]
t <- sub("_.*", "", TA); a <- sub(".*_", "", TA)
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set (see config.sh)"); v }
LAMBDA_OPT <- env("ARTP_LAMBDA"); N_OPT <- env("ARTP_N")
stopifnot(LAMBDA_OPT %in% c("1", "glgc"), N_OPT %in% c("median", "max"))
panel  <- c(eur = "eur", afr = "afr", eas = "eas", sas = "sas", his = "amr")[[a]]
suffix <- if (a %in% c("eur", "afr", "sas")) "HRC_1KGP3_others_ALL" else "1KGP3_ALL"
## the GLGC file sits in <trait>/ or, for the EUR runs of the multi-trait
## ablation (12_fetch_glgc.sh), in lipids/
gname  <- sprintf("%s_INV_%s_%s.meta.singlevar.results.gz", toupper(t), toupper(a), suffix)
glgc   <- file.path(env("GLGC_DATA"), c(t, "lipids"), gname)
glgc   <- if (any(file.exists(glgc))) glgc[file.exists(glgc)][1] else glgc[1]
genes_out <- file.path(env("GC_ROOT"), "inputs", "lipids", paste0(TA, ".genes.out"))
## the annotation the lipid MAGMA runs used (14_magma_trait_LSF.sh): EUR falls
## back to ldl_eur.genes.annot when g1000_eur.genes.annot is absent
annot     <- file.path(env("LIPID_ROOT"), paste0("g1000_", panel, ".genes.annot"))
if (panel == "eur" && !file.exists(annot)) annot <- file.path(env("LIPID_ROOT"), "ldl_eur.genes.annot")
bim       <- file.path(env("LIPID_ROOT"), paste0("g1000_", panel, ".bim"))
for (f in c(glgc, genes_out, annot, bim, env("REACTOME_SETS"))) if (!file.exists(f)) stop("missing ", f)
cat(sprintf("[%s] GLGC %s\n[%s] annotation %s\n", TA, glgc, TA, annot))
OUT <- env("ARTP_IN"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

## ---- pathways, exactly as CATFISH selected them ---------------------------------
g <- fread(genes_out, select = c("GENE", "P", "ZSTAT"))
tested <- as.character(g[is.finite(P) & is.finite(ZSTAT), GENE])
pw <- list()
for (ln in readLines(env("REACTOME_SETS"), warn = FALSE)) {
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]; if (length(f) < 2) next
  gs <- intersect(unique(f[-1][nzchar(f[-1])]), tested)
  if (length(gs) >= 10L && length(gs) <= 200L) pw[[f[1]]] <- gs
}
cat(sprintf("[%s] %d tested genes, %d pathways with 10-200 of them\n", TA, length(tested), length(pw)))

## ---- SNP -> gene from MAGMA's annotation ----------------------------------------
keep_genes <- unique(unlist(pw))
sg <- rbindlist(lapply(readLines(annot), function(ln) {
  if (startsWith(ln, "#")) return(NULL)
  f <- strsplit(ln, "\t", fixed = TRUE)[[1]]
  if (!(f[1] %in% keep_genes) || length(f) < 3) return(NULL)
  data.table(SNP = f[-(1:2)], Gene = f[1], Chr = as.integer(sub(":.*", "", f[2])))
}))
sg <- sg[!is.na(Chr) & Chr %in% 1:22]
cat(sprintf("[%s] %d SNP-gene pairs over %d genes\n", TA, nrow(sg), uniqueN(sg$Gene)))

## ---- summary statistics ---------------------------------------------------------
keep_file <- tempfile(fileext = ".txt"); writeLines(unique(sg$SNP), keep_file)
cmd <- sprintf("zcat %s | awk -F'\\t' -v OFS='\\t' 'NR==FNR {k[$1]=1; next} FNR==1 || ($1 in k)' %s -",
               shQuote(glgc), shQuote(keep_file))
ss <- fread(cmd = cmd, sep = "\t", select = c("rsID", "CHROM", "POS_b37", "REF", "ALT", "N",
                                              "EFFECT_SIZE", "SE", "pvalue_neg_log10", "pvalue_neg_log10_GC"),
            colClasses = list(numeric = c("EFFECT_SIZE", "SE", "pvalue_neg_log10", "pvalue_neg_log10_GC")))
unlink(keep_file)
inbim <- fread(bim, header = FALSE, select = 2L)[[1]]
ss <- ss[rsID %in% inbim & is.finite(EFFECT_SIZE) & is.finite(SE) & SE > 0]
ss <- ss[!duplicated(rsID)]
## GLGC's own lambda: chi-square of the raw p over that of their GC p (log scale,
## so no underflow), over SNPs with raw p > 1e-5
chi_lp <- function(nl) qchisq(-nl * log(10), 1, lower.tail = FALSE, log.p = TRUE)
lam_glgc <- with(ss[is.finite(pvalue_neg_log10) & is.finite(pvalue_neg_log10_GC) &
                    pvalue_neg_log10 > 0.301 & pvalue_neg_log10 < 5 & pvalue_neg_log10_GC > 0],
                 median(chi_lp(pvalue_neg_log10) / chi_lp(pvalue_neg_log10_GC)))
lambda  <- if (LAMBDA_OPT == "1") 1 else lam_glgc
nsample <- if (N_OPT == "median") median(ss$N) else max(ss$N)
summ <- ss[, .(SNP = rsID, RefAllele = REF, EffectAllele = ALT, BETA = EFFECT_SIZE, SE, P = pmax(10^(-pvalue_neg_log10), 1e-300),
               Chr = as.integer(CHROM), Pos = as.integer(POS_b37))]
sg <- sg[SNP %in% summ$SNP]
cat(sprintf("[%s] %d SNPs with statistics; lambda %s = %.4f (GLGC lambda %.4f); N %s = %d\n",
            TA, nrow(summ), LAMBDA_OPT, lambda, lam_glgc, N_OPT, as.integer(nsample)))

saveRDS(summ, file.path(OUT, paste0(TA, ".summary.rds")))
saveRDS(list(pathways = pw, snp_gene = sg), file.path(OUT, paste0(TA, ".pathways.rds")))
fwrite(data.table(set = TA, panel = panel, n_pathways = length(pw), n_snps = nrow(summ),
                  lambda_option = LAMBDA_OPT, lambda = lambda, lambda_glgc = lam_glgc,
                  n_option = N_OPT, nsamples = as.integer(nsample), glgc_file = glgc),
       file.path(OUT, paste0(TA, ".settings.tsv")), sep = "\t")
