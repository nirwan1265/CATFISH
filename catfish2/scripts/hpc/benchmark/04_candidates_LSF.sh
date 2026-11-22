#!/bin/bash
#BSUB -J "bench_genes[1-7]"
#BSUB -q sara
#BSUB -n 1
#BSUB -W 240
#BSUB -R "rusage[mem=16GB]"
#BSUB -o /share/maize/ntanduk/CATFISH/benchmark/logs/genes.%J.%I.out
#BSUB -e /share/maize/ntanduk/CATFISH/benchmark/logs/genes.%J.%I.err
############################################################################
## 04_candidates_LSF.sh -- the CATFISH gene score for the benchmark traits,
## with the paper's script (06_candidate_genes.R, copied to benchmark/scripts)
## and the lipid settings, only the input paths set by its environment
## variables:
##   MAGMA     magma_<trait>/<trait>_eur.genes.out
##   pathways  gc_real/out/<trait>_eur_reactome_5comp_westfall_young.csv (z_std),
##             Reactome sets
##   GWAS      <trait>/<trait>_eur.neglog10.txt via MAGMA's own .genes.annot
## Writes benchmark/candidates/06_candidate_genes_<trait>.csv
## Element i = TRAITS[i]. R env R_plot (needs data.table).
##   bsub -J "bench_genes[1-6]" < 04_candidates_LSF.sh    ; asthma: [7]
############################################################################
set -euo pipefail
TRAITS=(cad brca cd uc af sle asthma)
T=${TRAITS[$((LSB_JOBINDEX - 1))]}
ROOT=/share/maize/ntanduk/CATFISH
GC=$ROOT/gc_real
SCRIPT=$ROOT/benchmark/scripts/06_candidate_genes.R
OUTD=$ROOT/benchmark/candidates
export SPECIES=human TRAIT=$T
export MAGMA_GENES_OUT=$ROOT/magma_$T/${T}_eur.genes.out
export WY_CSV=$GC/out/${T}_eur_reactome_5comp_westfall_young.csv
export PATHWAY_SETS=$GC/inputs/c2.cp.reactome.v2026.2.Hs.entrez.sets.annot
export GENE_LOC=$ROOT/NCBI37.3.gene.loc
export MAGMA_GENES_ANNOT=$ROOT/ldl_eur.genes.annot
export GWAS_FILE=$ROOT/$T/${T}_eur.neglog10.txt
export CANDIDATE_OUT_DIR=$OUTD
export CATFISH_BASE=$ROOT/benchmark CATFISH_DATA=$ROOT
for f in "$SCRIPT" "$MAGMA_GENES_OUT" "$WY_CSV" "$PATHWAY_SETS" "$GENE_LOC" "$MAGMA_GENES_ANNOT" "$GWAS_FILE"; do
  [ -s "$f" ] || { echo "[$T] missing $f" >&2; exit 1; }
done
mkdir -p "$OUTD"
echo "[$T] $(date)"
Rscript "$SCRIPT"
echo "[$T] done $(date)"
