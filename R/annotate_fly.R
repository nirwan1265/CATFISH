################################################################################
### DGRP GWAS SNP -> gene annotation (window-based, NA-safe)
################################################################################

library(vroom)
library(dplyr)
library(GenomicRanges)
library(IRanges)
library(rtracklayer)
library(tibble)

################################################################################
### raw gwas
################################################################################
raw_gwas = vroom("/Users/nirwantandukar/Documents/Research/data/DGRP/Starvation_stress/raw_gwas/raw_GWAS_Starvation_stress_male_DGRP.csv")
head(raw_gwas)
colnames(raw_gwas)[4] <- "P.value"

################################################################################
### GFF3 file
################################################################################
gff  = "/Users/nirwantandukar/Documents/Research/data/Pathway/Drosophila_melanogaster.BDGP6.54.115.gff3"

ref_GRanges <- rtracklayer::import(gff)
genes_only  <- ref_GRanges[mcols(ref_GRanges)$type == "gene"]

# Ensure expected column names exist
# Your file looks like: CHR, Positions, P-Value, ...
if (!("CHR" %in% names(raw_gwas)))       stop("raw_gwas must have a CHR column")
if (!("Positions" %in% names(raw_gwas))) stop("raw_gwas must have a Positions column")
if (!("P.value" %in% names(raw_gwas)))   stop("raw_gwas must have a P.value column (rename P-Value to P.value)")

# SNP identifier: use SNP-Hash if present, else CHR:POS
snp_id <- if ("SNP-Hash" %in% names(raw_gwas)) raw_gwas[["SNP-Hash"]] else paste0(raw_gwas$CHR, ":", raw_gwas$Positions)

# Filter criteria (change as needed)
p_cut <- 0.05   # e.g., 1e-5, 1e-6, 0.05, etc.
keep <- !is.na(raw_gwas$P.value) & raw_gwas$P.value <= p_cut

# Build SNP GRanges (Drosophila chr names already match GFF seqnames usually)
snps <- GRanges(
  seqnames = Rle(as.character(raw_gwas$CHR)),
  ranges   = IRanges(start = raw_gwas$Positions, end = raw_gwas$Positions)
)
mcols(snps)$SNP     <- as.character(snp_id)
mcols(snps)$P.value <- raw_gwas$P.value

snps_f <- snps[keep]
if (length(snps_f) == 0) stop("No SNPs passed filtering (p_cut too strict?)")

# Window around SNP (±12.5kb => width 25001 like your example)
win_width <- 25001
snps_w <- resize(snps_f, width = win_width, fix = "center")

# Optional: restrict to standard fly chroms to avoid weird contigs
std_chr <- c("2L","2R","3L","3R","4","X","Y")
snps_w  <- snps_w[as.character(seqnames(snps_w)) %in% std_chr]
genes_only2 <- genes_only[as.character(seqnames(genes_only)) %in% std_chr]

# Overlaps
ol <- findOverlaps(snps_w, genes_only2, ignore.strand = TRUE)
if (length(ol) == 0) stop("No overlaps found. Check chromosome naming between GWAS and GFF.")

hs <- queryHits(ol)   # SNP window hits
hg <- subjectHits(ol) # gene hits

# Relation of SNP *position* to gene body
# (Use center of window = original SNP position)
snp_pos <- start(snps_f)[keep][hs]  # careful: snps_f is already filtered; indices align with snps_w
# safer: store original SNP position as metadata before resizing:
mcols(snps_w)$SNP_pos <- start(snps_f)

snp_pos <- mcols(snps_w)$SNP_pos[hs]
gene_start <- start(genes_only2)[hg]
gene_end   <- end(genes_only2)[hg]

ann <- tibble(
  GeneID   = mcols(genes_only2)$ID[hg],
  GeneName = if ("Name" %in% names(mcols(genes_only2))) mcols(genes_only2)$Name[hg] else NA_character_,
  CHR      = as.character(seqnames(genes_only2))[hg],
  SNP      = mcols(snps_w)$SNP[hs],
  Pos      = snp_pos,
  P.value  = mcols(snps_w)$P.value[hs],
  Relation = ifelse(
    snp_pos >= gene_start & snp_pos <= gene_end, "within",
    ifelse(snp_pos < gene_start, "upstream", "downstream")
  ),
  DistToGene = ifelse(
    snp_pos < gene_start, gene_start - snp_pos,
    ifelse(snp_pos > gene_end, snp_pos - gene_end, 0)
  )
) %>%
  arrange(P.value, CHR, Pos)

ann

# Save
write.csv(ann,"Individual_GWAS_male_starvation.csv", row.names = F)
