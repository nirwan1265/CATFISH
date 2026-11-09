#!/usr/bin/env Rscript
## =============================================================================
## 34_artp_collect.R -- gather each lipid set's sARTP chunks into
## ARTP_OUT/<ta>.csv and report pathways missing or failed.
##   Rscript 34_artp_collect.R          (or bsub < 34_artp_collect_LSF.sh)
## =============================================================================
suppressPackageStartupMessages(library(data.table))
IN <- Sys.getenv("ARTP_IN"); OUT <- Sys.getenv("ARTP_OUT")
for (TA in strsplit(Sys.getenv("LIPIDS"), " ")[[1]]) {
  pf <- file.path(IN, paste0(TA, ".pathways.rds")); if (!file.exists(pf)) next
  want <- names(readRDS(pf)$pathways)
  fs <- list.files(file.path(OUT, TA), pattern = "^chunk[0-9]+\\.csv$", full.names = TRUE)
  r <- if (length(fs)) unique(rbindlist(lapply(fs, fread), fill = TRUE)[order(!is.finite(p_artp))], by = "pathway_id") else data.table()
  fwrite(r, file.path(OUT, paste0(TA, ".csv")))
  ok <- if (nrow(r)) sum(is.finite(r$p_artp)) else 0L
  cat(sprintf("%-8s %4d / %4d pathways with p, %d failed, %d missing; Bonferroni-significant %d\n",
              TA, ok, length(want), if (nrow(r)) sum(!is.finite(r$p_artp)) else 0L,
              length(setdiff(want, r$pathway_id)),
              if (nrow(r)) sum(r$p_artp < 0.05 / length(want), na.rm = TRUE) else 0L))
}
