#!/usr/bin/env Rscript
## 08a_write_set_annot.R -- the SorghumBicolorCyc pathways as a MAGMA
## --set-annot file, read exactly as 05_wy_perm.R reads them (Pathway-id,
## Gene-name; "unknown" and empty genes dropped), so that MAGMA's gene-set
## analysis sees the same pathways as CATFISH. Size limits are applied later
## (09), on the genes present, as in CATFISH.
##   Rscript 08a_write_set_annot.R          -> INPUT_DIR/sorghumbicolorcyc.sets.annot
suppressPackageStartupMessages(library(CATFISH))
getcfg <- function(k) { v <- Sys.getenv(k); if (v == "") stop("config missing: ", k); v }
pw_raw <- utils::read.delim(getcfg("PATHWAY_FILE"), stringsAsFactors = FALSE, check.names = FALSE)
pw <- unique(data.frame(pathway_id = pw_raw[["Pathway-id"]], gene_id = pw_raw[["Gene-name"]],
                        stringsAsFactors = FALSE))
pw <- pw[nzchar(pw$gene_id) & pw$gene_id != "unknown", ]
out <- file.path(getcfg("INPUT_DIR"), "sorghumbicolorcyc.sets.annot")
CATFISH:::catfish_write_set_annot(pw, out)
cat(sprintf("%d pathways -> %s\n", length(unique(pw$pathway_id)), out))
