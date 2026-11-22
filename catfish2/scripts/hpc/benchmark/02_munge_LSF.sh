#!/bin/bash
#BSUB -J "bench_munge[1-7]"
#BSUB -n 1
#BSUB -W 180
##BSUB -q sara
#BSUB -R "rusage[mem=8GB]"
#BSUB -o /share/maize/ntanduk/CATFISH/benchmark/logs/munge.%J.%I.out
#BSUB -e /share/maize/ntanduk/CATFISH/benchmark/logs/munge.%J.%I.err
############################################################################
## 02_munge_LSF.sh -- GWAS Catalog harmonised summary statistics -> the two
## files the lipid pipeline reads, one array element per traits.tsv line
## (element i = data line i). Same outputs and filters as 13_munge_trait.sh:
##   /share/maize/ntanduk/CATFISH/<trait>/<trait>_eur.magma.pval    SNP P N
##                     (p floored at 1e-300 for MAGMA)
##   /share/maize/ntanduk/CATFISH/<trait>/<trait>_eur.neglog10.txt  SNP NEGLOG10P
##                     (unfloored, for the candidate-gene GWAS layer)
## Variants: rsIDs (hm_rsid) present in the 1000G EUR panel used for MAGMA,
## first occurrence kept. N: the file's per-variant sample size column if it
## has one (Neff, n, N, sample_size), otherwise n_total from traits.tsv.
## Columns are found by header name, so the GWAS Catalog layouts all work.
## Then MAGMA exactly as the lipids: TRAIT=<trait> bsub < ../../14_magma_trait_LSF.sh
##   bsub < 02_munge_LSF.sh                (all 7)
##   bsub -J "bench_munge[2]" < 02_munge_LSF.sh   (one trait: line 2 = cad)
############################################################################
set -euo pipefail
SCRIPTS=/share/maize/ntanduk/CATFISH/benchmark/scripts
SUMS=/share/maize/ntanduk/CATFISH/benchmark/sumstats
ROOT=/share/maize/ntanduk/CATFISH
BIM=$ROOT/g1000_eur.bim

LINE=$(awk -F'\t' -v i="$LSB_JOBINDEX" 'NR == i + 1' "$SCRIPTS/traits.tsv")
[ -n "$LINE" ] || { echo "no line $LSB_JOBINDEX in traits.tsv" >&2; exit 1; }
T=$(echo "$LINE" | cut -f1); NTOT=$(echo "$LINE" | cut -f5)
IN=$(ls "$SUMS/$T"/*.h.tsv.gz 2>/dev/null | head -1 || true)
[ -n "$IN" ] && [ -f "$IN" ] || { echo "[$T] no .h.tsv.gz in $SUMS/$T" >&2; exit 1; }
[ -f "$BIM" ] || { echo "missing $BIM" >&2; exit 1; }
mkdir -p "$ROOT/$T"
PVAL=$ROOT/$T/${T}_eur.magma.pval
NLOG=$ROOT/$T/${T}_eur.neglog10.txt
echo "[$T] $(date) in=$IN  n_total=$NTOT"

zcat "$IN" | awk -F'\t' -v OFS='\t' -v BIM="$BIM" -v PVAL="$PVAL" -v NLOG="$NLOG" -v NTOT="$NTOT" '
  function nlog10(s,   p, m, e, k) {          # -log10 of a p-value string, unfloored
    p = s + 0
    if (p > 0) return -log(p) / log(10)
    k = split(tolower(s), m, "e")             # underflow (e.g. 1e-400): mantissa and exponent
    if (k == 2 && (m[1] + 0) > 0) return -(log(m[1] + 0) / log(10) + (m[2] + 0))
    return -1
  }
  BEGIN {
    while ((getline line < BIM) > 0) { split(line, b, "\t"); ref[b[2]] = 1; nref++ }
    close(BIM)
    printf("reference SNPs loaded: %d\n", nref) > "/dev/stderr"
    print "SNP", "P", "N"    > PVAL
    print "SNP", "NEGLOG10P" > NLOG
  }
  NR == 1 {
    for (i = 1; i <= NF; i++) h[$i] = i
    cs = ("hm_rsid" in h) ? h["hm_rsid"] : (("variant_id" in h) ? h["variant_id"] : 0)
    cp = ("p_value" in h) ? h["p_value"] : 0
    cn = 0; split("Neff n N sample_size n_total", nn, " ")
    for (j = 1; j <= 5; j++) if (nn[j] in h) { cn = h[nn[j]]; break }
    if (!cs || !cp) { print "no hm_rsid/variant_id or p_value column" > "/dev/stderr"; exit 1 }
    printf("columns: SNP=%s P=%s N=%s\n", $cs, $cp, cn ? $cn : "n_total " NTOT) > "/dev/stderr"
    next
  }
  $cs !~ /^rs/    { nonrs++;  next }
  !($cs in ref)   { notref++; next }
  ($cs in seen)   { dup++;    next }
  {
    nl = nlog10($cp)
    if (nl < 0 || nl != nl) { bad++; next }
    n = cn ? $cn + 0 : NTOT + 0
    if (n <= 0) { badn++; next }
    seen[$cs] = 1
    if (nl > maxnl) maxnl = nl
    if (nl > 300) { p = 1e-300; floored++ } else p = 10 ^ (-nl)
    printf "%s\t%.6g\t%d\n", $cs, p, n  > PVAL
    printf "%s\t%.6f\n",     $cs, nl    > NLOG
    kept++
  }
  END {
    printf("kept=%d not_in_panel=%d non_rs=%d duplicate=%d bad_p=%d bad_n=%d floored=%d max_neglog10=%.1f\n",
           kept, notref, nonrs, dup, bad, badn, floored, maxnl) > "/dev/stderr"
  }
'
echo "[$T] wrote $PVAL ($(( $(wc -l < "$PVAL") - 1 )) variants)"
echo "[$T] wrote $NLOG"
echo "[$T] done $(date)"
