#!/usr/bin/env Rscript
## =============================================================================
## 03_write_pval.R <rep_id> <chr|all> <magma|pascalx> <out file>
## One replicate's SNP p-values (two-sided, from its z) in the format the lipid
## runs used: "SNP P N" for MAGMA (N = SIM_N), "SNP P" for PascalX. Used by
## steps 03 and 07 to write a temporary input file.
## =============================================================================
suppressPackageStartupMessages(library(data.table))
env <- function(k) { v <- Sys.getenv(k); if (!nzchar(v)) stop(k, " is not set"); v }
a <- commandArgs(trailingOnly = TRUE)
r <- as.integer(a[1]); CH <- a[2]; FMT <- a[3]; OUTF <- a[4]
RUN <- env("RUN")
d <- fread(file.path(RUN, "design.tsv"))[rep_id == r]
if (nrow(d) != 1) stop("replicate ", r, " is not in the design")
col <- sprintf("rep_%04d", r)
z <- fread(file.path(RUN, "z", sprintf("z_cell%02d.tsv.gz", d$cell)), select = c("SNP", "CHR", col))
if (CH != "all") z <- z[CHR == as.integer(CH)]
p <- 2 * stats::pnorm(-abs(z[[col]]))
o <- if (FMT == "magma") data.table(SNP = z$SNP, P = p, N = as.integer(env("SIM_N"))) else
     data.table(SNP = z$SNP, P = p)
fwrite(o, OUTF, sep = "\t")
