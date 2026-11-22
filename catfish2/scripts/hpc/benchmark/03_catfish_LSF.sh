#!/bin/bash
#BSUB -J "bench_catfish[1-7]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 720
#BSUB -o /share/maize/ntanduk/CATFISH/benchmark/logs/catfish.%J.%I.out
#BSUB -e /share/maize/ntanduk/CATFISH/benchmark/logs/catfish.%J.%I.err
############################################################################
## 03_catfish_LSF.sh -- CATFISH Westfall-Young pathway analysis for the
## benchmark traits, exactly as for the lipids (same script, same settings):
##   1. gene-gene correlations from the trait's merged .genes.raw, as
##      scripts/human_ldl/17_prep_trait_catfish.R (Entrez IDs, all pairs)
##      -> gc_real/inputs/lipids/<trait>_eur.cor_pairs.txt
##   2. the merged .genes.out copied to gc_real/inputs/lipids/<trait>_eur.genes.out
##   3. gc_real/wy_gc_real.R <trait>_eur  (Reactome, 10-200 genes, B = 10000,
##      5 components, tau {0.01, 0.05, 0.5}, seed 20260706, gene-level GC)
##      -> gc_real/out/<trait>_eur_reactome_5comp_westfall_young.csv
## Element i = TRAITS[i]. Submit with the R env that has CATFISH (R_plot).
##   bsub < 03_catfish_LSF.sh                     (cad brca cd uc af sle; asthma = 7)
##   bsub -J "bench_catfish[7]" < 03_catfish_LSF.sh
############################################################################
set -euo pipefail
TRAITS=(cad brca cd uc af sle asthma)
T=${TRAITS[$((LSB_JOBINDEX - 1))]}
ROOT=/share/maize/ntanduk/CATFISH
GC=$ROOT/gc_real
RAW=$ROOT/magma_$T/${T}_eur.genes.raw
GOUT=$ROOT/magma_$T/${T}_eur.genes.out
IN=$GC/inputs/lipids
for f in "$RAW" "$GOUT" "$GC/wy_gc_real.R" "$GC/inputs/c2.cp.reactome.v2026.2.Hs.entrez.sets.annot"; do
  [ -s "$f" ] || { echo "[$T] missing $f" >&2; exit 1; }
done
mkdir -p "$IN"
echo "[$T] $(date) correlations from $RAW"

Rscript -e '
suppressPackageStartupMessages(library(CATFISH))
a <- commandArgs(trailingOnly = TRUE)
magma_genesraw_to_cor_pairs_banded(genes_raw_file = a[1], out_pairs_file = a[2],
  gene_regex = "^[0-9]", keep_abs_r_ge = 0, overwrite = TRUE, verbose = FALSE)
g  <- utils::read.table(a[3], header = TRUE, stringsAsFactors = FALSE)
cp <- utils::read.delim(a[2], header = TRUE, stringsAsFactors = FALSE)
chr <- stats::setNames(as.character(g$CHR), as.character(g$GENE))
cov <- unique(stats::na.omit(chr[as.character(cp[[1]])]))
cat(sprintf("%d gene pairs; %d genes across %d chromosomes covered\n", nrow(cp), nrow(g), length(cov)))
if (length(cov) < 22) stop("correlations cover fewer than 22 chromosomes")
' "$RAW" "$IN/${T}_eur.cor_pairs.txt" "$GOUT"

cp "$GOUT" "$IN/${T}_eur.genes.out"
cd "$GC"
GC_ROOT="$GC" Rscript "$GC/wy_gc_real.R" "${T}_eur"
echo "[$T] done $(date)"
