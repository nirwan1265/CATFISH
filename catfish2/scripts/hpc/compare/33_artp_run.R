#!/usr/bin/env Rscript
## =============================================================================
## 33_artp_run.R -- sARTP (ARTP2::sARTP, method = 3, the ARTP test) for one
## chunk of one lipid set's pathways. Each pathway is one sARTP call with:
##   summary.files  that pathway's SNPs from <ta>.summary.rds (written to a
##                  small file in the job's own work directory)
##   pathway        SNP / Gene / Chr from MAGMA's SNP-to-gene map
##   reference      the set's 1000 Genomes panel, one row per chromosome
##   family         gaussian, lambda and nsamples from <ta>.settings.tsv
##   options        nperm, maf, min.marg.p from config.sh; everything else at
##                  sARTP's defaults; nthread = the cores LSF gave the job
## Pathways already in the chunk's output with a p-value are skipped, so a killed
## job resumes and a rerun retries only the ones that failed.
##   Rscript 33_artp_run.R <set> <chunk>      chunk in 1..ARTP_NCHUNK
## =============================================================================
suppressPackageStartupMessages(library(data.table))
.libPaths(c(Sys.getenv("ARTP_RLIB"), .libPaths()))
suppressPackageStartupMessages(library(ARTP2))
args <- commandArgs(trailingOnly = TRUE); TA <- args[1]; K <- as.integer(args[2])
IN <- Sys.getenv("ARTP_IN"); NCH <- as.integer(Sys.getenv("ARTP_NCHUNK"))
st <- fread(file.path(IN, paste0(TA, ".settings.tsv")))
summ <- readRDS(file.path(IN, paste0(TA, ".summary.rds")))
pp   <- readRDS(file.path(IN, paste0(TA, ".pathways.rds")))
ids  <- sort(names(pp$pathways)); ids <- ids[seq(K, length(ids), by = NCH)]
setkey(summ, SNP); sg <- pp$snp_gene; setkey(sg, Gene)

OUTD <- file.path(Sys.getenv("ARTP_OUT"), TA); dir.create(OUTD, recursive = TRUE, showWarnings = FALSE)
outf <- file.path(OUTD, sprintf("chunk%02d.csv", K))
done <- if (file.exists(outf)) fread(outf)[is.finite(p_artp), pathway_id] else character(0)   # failed ones are retried
ids  <- setdiff(ids, done)
WORK <- file.path(OUTD, sprintf("work_chunk%02d", K)); dir.create(WORK, showWarnings = FALSE)
setwd(WORK)                                     # sARTP's own temporary files stay per job

ref_all <- data.frame(bed = file.path(Sys.getenv("ARTP_REF"), sprintf("g1000_%s.chr%d.bed", st$panel, 1:22)),
                      bim = file.path(Sys.getenv("ARTP_REF"), sprintf("g1000_%s.chr%d.bim", st$panel, 1:22)),
                      fam = file.path(Sys.getenv("ARTP_REF"), sprintf("g1000_%s.chr%d.fam", st$panel, 1:22)),
                      stringsAsFactors = FALSE)
stopifnot(all(file.exists(ref_all$bed)))
nthr <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "1"))
opts <- list(method = 3, nperm = as.integer(as.numeric(Sys.getenv("ARTP_NPERM"))),
             maf = as.numeric(Sys.getenv("ARTP_MAF")),
             min.marg.p = as.numeric(Sys.getenv("ARTP_MIN_MARG_P")),
             nthread = nthr)
cat(sprintf("[%s chunk %d] %d pathways to run (%d done), lambda %.4f, N %d, nperm %d, maf %g, min.marg.p %g, %d threads\n",
            TA, K, length(ids), length(done), st$lambda, st$nsamples, opts$nperm, opts$maf, opts$min.marg.p, nthr))

for (id in ids) {
  t0 <- Sys.time()
  pwd <- unique(sg[J(pp$pathways[[id]]), nomatch = 0L][, .(SNP, Gene, Chr)])
  chrs <- sort(unique(pwd$Chr))
  sf <- file.path(WORK, "summary.txt")
  fwrite(summ[J(unique(pwd$SNP)), nomatch = 0L], sf, sep = "\t")
  r <- tryCatch(sARTP(summary.files = sf, pathway = as.data.frame(pwd), family = "gaussian",
                      reference = ref_all[chrs, , drop = FALSE], lambda = st$lambda,
                      nsamples = list(st$nsamples), options = opts),
                error = function(e) e)
  row <- if (inherits(r, "error")) {
    data.table(pathway_id = id, n_genes_in = length(pp$pathways[[id]]), n_snps_in = uniqueN(pwd$SNP),
               n_genes_used = NA_integer_, p_artp = NA_real_, accurate = NA,
               error = gsub("[\r\n,]", " ", conditionMessage(r)))
  } else {
    data.table(pathway_id = id, n_genes_in = length(pp$pathways[[id]]), n_snps_in = uniqueN(pwd$SNP),
               n_genes_used = if (is.data.frame(r$gene.pvalue)) nrow(r$gene.pvalue) else NA_integer_,
               p_artp = r$pathway.pvalue, accurate = isTRUE(r$accurate), error = "")
  }
  row[, seconds := round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)]
  fwrite(row, outf, append = file.exists(outf))
  cat(sprintf("  %s  p = %s  (%.0f s)\n", id, format(row$p_artp, digits = 3), row$seconds))
}
cat(sprintf("[%s chunk %d] done -> %s\n", TA, K, outf))
