#!/usr/bin/env Rscript
# =============================================================================
# 10_spike_setup.R  --  Pick the causal genes/SNPs for the staged spike-in.
#
# Takes the target pathway, finds its genes, looks them up in the EXISTING
# per-chromosome MAGMA .genes.annot to get their SNPs, and picks ONE
# representative causal SNP per causal gene. Writes:
#   ${SPIKE_WORK}/causal_snps.txt   (plain SNP id list, for PLINK --extract)
#   ${SPIKE_WORK}/causal_map.tsv    (SNP <-> GENE map, for the simulator)
# =============================================================================
suppressMessages({
  library(data.table)
})

getcfg <- function(k) {
  v <- Sys.getenv(k, "")
  if (!nzchar(v)) stop(paste("config missing:", k), call. = FALSE)
  v
}

SPECIES <- getcfg("SPECIES")
SPIKE_PATHWAY <- getcfg("SPIKE_PATHWAY")
ANNOT_TMPL <- getcfg("ANNOT_PATH_TEMPLATE")
CHROMS <- strsplit(trimws(getcfg("CHROMS")), "\\s+")[[1]]
SPIKE_WORK <- getcfg("SPIKE_WORK")
CATFISH_REPO <- getcfg("CATFISH_REPO")
PATHWAY_FILE <- getcfg("PATHWAY_FILE")
SPIKE_GENE_COL <- Sys.getenv("SPIKE_GENE_COL", "")
set.seed(as.integer(getcfg("SPIKE_SEED")))

pw_raw <- fread(PATHWAY_FILE, check.names = FALSE)
pw_sub <- pw_raw[`Pathway-id` == SPIKE_PATHWAY]
if (!nrow(pw_sub)) {
  stop(sprintf("Pathway '%s' not found / has no genes. Set SPIKE_PATHWAY to a real pathway_id.", SPIKE_PATHWAY), call. = FALSE)
}
cat(sprintf("[spike] pathway %s found with %d rows in %s\n", SPIKE_PATHWAY, nrow(pw_sub), basename(PATHWAY_FILE)))

gene2snps <- list()
for (chr in CHROMS) {
  af <- gsub("@CHR@", chr, ANNOT_TMPL, fixed = TRUE)
  if (!file.exists(af)) {
    warning(paste("annot missing:", af))
    next
  }
  ln <- readLines(af, warn = FALSE)
  ln <- ln[!startsWith(ln, "#") & nzchar(ln)]
  for (l in ln) {
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    g <- f[1]
    snps <- if (length(f) >= 3) f[3:length(f)] else character(0)
    if (length(snps)) gene2snps[[g]] <- snps
  }
}
cat(sprintf("[spike] annotation covers %d genes total\n", length(gene2snps)))

annot_keys <- toupper(names(gene2snps))
pathway_candidates <- intersect(c("Gene-name", "Gene-id"), names(pw_sub))
if (nzchar(SPIKE_GENE_COL)) {
  if (!SPIKE_GENE_COL %in% names(pw_sub)) {
    stop(sprintf("SPIKE_GENE_COL='%s' is not a column in %s", SPIKE_GENE_COL, basename(PATHWAY_FILE)), call. = FALSE)
  }
  pathway_candidates <- c(SPIKE_GENE_COL, setdiff(pathway_candidates, SPIKE_GENE_COL))
}
if (!length(pathway_candidates)) {
  stop("Pathway file lacks both 'Gene-name' and 'Gene-id' columns.", call. = FALSE)
}

overlap_tbl <- rbindlist(lapply(pathway_candidates, function(col) {
  vals <- unique(as.character(pw_sub[[col]]))
  vals <- vals[nzchar(vals) & !is.na(vals)]
  data.table(
    gene_col = col,
    n_pathway_genes = length(vals),
    n_overlap = sum(toupper(vals) %in% annot_keys)
  )
}))
setorder(overlap_tbl, -n_overlap, -n_pathway_genes)
print(overlap_tbl)

best_col <- overlap_tbl$gene_col[[1]]
genes_in_pw <- unique(as.character(pw_sub[[best_col]]))
genes_in_pw <- genes_in_pw[nzchar(genes_in_pw) & !is.na(genes_in_pw)]
cat(sprintf("[spike] using %s for pathway gene IDs (%d genes; %d overlap with annot)\n",
            best_col, overlap_tbl$n_pathway_genes[[1]], overlap_tbl$n_overlap[[1]]))

gmatch <- intersect(toupper(genes_in_pw), annot_keys)
name_lookup <- setNames(names(gene2snps), annot_keys)
causal_genes <- unname(name_lookup[gmatch])
if (!length(causal_genes)) {
  ex_path <- paste(head(genes_in_pw, 5), collapse = ", ")
  ex_annot <- paste(head(names(gene2snps), 5), collapse = ", ")
  stop(
    "None of the pathway genes have annotated SNPs.\n",
    "Pathway examples: ", ex_path, "\n",
    "Annot examples: ", ex_annot, "\n",
    "Try setting SPIKE_GENE_COL to 'Gene-name' or 'Gene-id' in config.spikein.sh.",
    call. = FALSE
  )
}
cat(sprintf("[spike] %d of %d pathway genes have SNPs and will be causal\n",
            length(causal_genes), length(genes_in_pw)))

map <- rbindlist(lapply(causal_genes, function(g) {
  snps <- gene2snps[[g]]
  data.table(GENE = g, SNP = snps[sample.int(length(snps), 1)])
}))

dir.create(SPIKE_WORK, recursive = TRUE, showWarnings = FALSE)
fwrite(map, file.path(SPIKE_WORK, "causal_map.tsv"), sep = "\t")
writeLines(map$SNP, file.path(SPIKE_WORK, "causal_snps.txt"))
cat(sprintf("[spike] wrote %d causal SNPs -> %s\n",
            nrow(map), file.path(SPIKE_WORK, "causal_snps.txt")))
