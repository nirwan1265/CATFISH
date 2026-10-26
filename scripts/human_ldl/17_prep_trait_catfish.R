#!/usr/bin/env Rscript
## =============================================================================
## 17_prep_trait_catfish.R -- build the gene-gene correlation pairs for one
## GLGC lipid trait from its merged MAGMA .genes.raw.
##
## This is the step between the MAGMA merge and 18_westfall_young_trait.R.
## It must run on the MERGED .genes.raw, not on the per-chromosome batches:
## the merged file is what carries correlations for every chromosome, and a
## batch file carries only its own.
##
##   TRAIT=hdl Rscript 17_prep_trait_catfish.R
## =============================================================================
suppressPackageStartupMessages({ library(CATFISH) })

TRAIT <- tolower(Sys.getenv("TRAIT", "hdl"))
BASE  <- path.expand(Sys.getenv("CATFISH_BASE", "~/Documents/Github/MAGCAT"))
MDIR  <- Sys.getenv("MAGMA_DIR",
                    file.path(BASE, "results/human_LDL", paste0("magma_", TRAIT)))
RAW   <- file.path(MDIR, sprintf("%s_eur.genes.raw", TRAIT))
OUT_D <- file.path(BASE, "results/human_LDL/catfish")
OUT   <- file.path(OUT_D, sprintf("%s_eur.cor_pairs.txt", TRAIT))
dir.create(OUT_D, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(RAW)) {
  batches <- list.files(MDIR, pattern = "batch[0-9]+_chr\\.genes\\.raw$")
  stop(sprintf(paste0("merged file not found: %s\n",
       "  %d per-chromosome batch files are present, but the merge step has\n",
       "  not been run (or its output was not copied back). On the cluster:\n",
       "    TRAIT=%s bsub < 15_magma_merge_trait_LSF.sh"),
       RAW, length(batches), TRAIT), call. = FALSE)
}

cat(sprintf("trait=%s\nreading %s\n", TRAIT, RAW))
## argument names and defaults as in R/correlation_genes_wrappers.R --
## gene_regex defaults to "^FBgn" (FlyBase) and silently matches nothing for
## human Entrez IDs, so it must be set explicitly.
magma_genesraw_to_cor_pairs_banded(
  genes_raw_file = RAW,
  out_pairs_file = OUT,
  gene_regex     = "^[0-9]",   # human Entrez IDs
  keep_abs_r_ge  = 0,
  overwrite      = TRUE,
  verbose        = TRUE
)
cat(sprintf("wrote %s\n", OUT))

## sanity check: how many chromosomes are actually represented?
go <- file.path(MDIR, sprintf("%s_eur.genes.out", TRAIT))
if (file.exists(go)) {
  g <- utils::read.table(go, header = TRUE, stringsAsFactors = FALSE)
  cp <- utils::read.delim(OUT, header = TRUE, stringsAsFactors = FALSE)
  names(cp)[1:2] <- c("g1", "g2")
  chr <- stats::setNames(as.character(g$CHR), as.character(g$GENE))
  covered <- unique(chr[as.character(cp$g1)])
  covered <- covered[!is.na(covered)]
  cat(sprintf("genes analysed: %d across %d chromosomes\n",
              nrow(g), length(unique(g$CHR))))
  cat(sprintf("correlation pairs: %d covering %d chromosomes\n",
              nrow(cp), length(covered)))
  if (length(covered) < length(unique(g$CHR))) {
    warning(sprintf(paste0("correlations cover %d of %d chromosomes. ",
      "Genes on the rest will be treated as INDEPENDENT and the LD-aware ",
      "null will be incomplete."), length(covered), length(unique(g$CHR))),
      call. = FALSE)
  }
}
