# simulation_check_catfish — CATFISH vs MAGMA vs PascalX with known truth

Pipeline written Nov 2026; nothing has been run yet. Pilot first, then the full run.

## Question
On real lipid GWAS all methods put plasma lipoprotein first, and the competitive
calls are of similar size (LDL/HDL EUR: CATFISH 9/14, MAGMA competitive 5/12,
PascalX 13/17). Real data cannot say which method is more powerful or which keeps
false positives at the stated rate. This simulation puts signal in known genes of
known pathways, on real LD, and scores every method on the same GWAS statistics.

## Settings
- LD: the 1000 Genomes EUR panel of the lipid runs (g1000_eur, 503 people)
- Chromosomes 1-10; HapMap3 SNPs (LDSC w_hm3.snplist); MAF >= 0.01 in the panel
- Genes: MAGMA's EUR annotation (+/- 10 kb, NCBI37.3, Entrez IDs) restricted to
  the panel's SNPs. Pathways: Reactome restricted to those genes, 10-200 genes.
- 50 replicates per scenario x target size (pilot: 5)
- CATFISH: catfish_wy_pathways(), B = 2,000 for both nulls, five components,
  tau {0.01, 0.05, 0.5}, seed 20260706, genomic control on (package default),
  competitive test on (genes in CHR/START order)
- MAGMA: default gene model (as the lipid runs), gene-set analysis competitive
  and --model self-contained. PascalX: +/- 10 kb, MAF 0.01 (as compare/21).
- Calls at family-wise 0.05: CATFISH Westfall-Young; MAGMA and PascalX
  Bonferroni over the pathways tested. sARTP is not included.

## How the GWAS statistics are made
Noise: for standardised genotypes X (503 x M) and u = e / sqrt(503),
e ~ N(0, I), z = X'u has exactly the panel's LD as its correlation. Step 01
draws u; step 02 gets X'u for every panel SNP from plink2 --variant-score,
standardises it, and checks it against the exact z of the causal SNPs (up to
the sign set by which allele is counted; the job stops otherwise).
Signal: step 02 adds each causal SNP's expected z (lambda, random sign) to the
SNPs within LD_WINDOW_KB (1 Mb), weighted by their panel correlation. Local LD
only: with 503 people, unlinked SNPs show chance correlations of about
+/- 0.045, and spreading the signal through them inflated every SNP (genome-wide
SD of z 2.1-2.8 in the first pilot attempt).

## Scenarios
| code | truth |
|---|---|
| S0 | no causal SNPs |
| S1 | background only: 2% of genes, one causal SNP each, expected z 3 (random sign) |
| S2 SSS | background + target pathway: 15% of its genes at delta 3.00 |
| S3 MDS | background + 60% of target genes at delta 1.10 |
| S4 WDE | background + 90% of target genes at delta 0.55 |
| S5 MSA | background + 10% at 1.80 and 60% at 1.10 |
| S6 SGP | background + 5% (at least one gene) at 4.00 |

Layers as in 17_archetype_power_wy.R. Each causal gene gets one causal SNP with
expected z = delta x EFFECT_SCALE (random sign). Targets: one fixed non-lipid
Reactome pathway per size class (12-18, 35-45, 90-110 genes on chr 1-10).
Background genes exclude the target's genes.

## What is measured (08_aggregate.R)
- power: the target pathway is called (S2-S6)
- off-target calls: pathways called with no gene within 1 Mb of a target causal
  SNP (in S0 and S1, every call). Read as an error where it is one:
  S0 all methods; S1 and S2-S6 the competitive tests and "both". The
  self-contained tests are expected to call pathways carrying background genes.

## Steps (separate LSF submissions; no memory request)
No script activates a conda env: submit the R steps with the R env active
and 00b / 07b with the PascalX env active.
| step | script | what |
|---|---|---|
| 00 | 00_prepare_panel_LSF.sh (+ 00_prepare_sets.R) | panel, genotype counts, annotation, pathways |
| 00b | 00b_pascalx_ref_LSF.sh | PascalX import of the panel |
| 01 | 01_make_scenarios_LSF.sh (+ .R) | design, targets, causal SNPs, u |
| 02 | 02_simulate_z_LSF.sh (+ .R) | z for every replicate (one element per design cell) |
| 03 | 03_magma_gene_LSF.sh (+ 03_write_pval.R) | MAGMA gene analysis, one element per replicate x chromosome |
| 03b | 03b_magma_merge_LSF.sh | merge the 10 chromosomes of each replicate |
| 04 | 04_corpairs_once_LSF.sh (+ .R) | gene-gene correlations, once (LD only) |
| 05 | 05_catfish_LSF.sh (+ .R) | CATFISH, one element per replicate |
| 06 | 06_magma_gsa_LSF.sh | MAGMA gene-set, one element per replicate |
| 07a | 07a_pascalx_pval_LSF.sh | PascalX input p-values (R env), one element per replicate |
| 07b | 07b_pascalx_LSF.sh | PascalX (PascalX env), one element per replicate |
| 08 | 08_aggregate_LSF.sh (+ .R) | power and off-target rates |

The simulated statistics (02) and MAGMA (03) are separate submissions. Steps 03
and 07 write each replicate's SNP p-values to a temporary file and remove it
when the job ends; step 02 removes plink2's raw .vscore once z is written.

## Before the pilot
1. Copy this folder to /share/maize/ntanduk/CATFISH/simulation_check_catfish
   and w_hm3.snplist to its inputs/.
2. CATFISH with the competitive test: either reinstall the package, or keep
   DEV_WY (config.sh) pointing at the new westfall_young.R, as for the EUR test.

## Pilot (MODE=pilot, the default): S0 and S3 (medium target) at scales 2, 4, 6
4 cells x 5 replicates = 20 replicates.
```
cd /share/maize/ntanduk/CATFISH/simulation_check_catfish
bsub < 00_prepare_panel_LSF.sh                      # then
bsub < 00b_pascalx_ref_LSF.sh                       # (PascalX env)
bsub < 01_make_scenarios_LSF.sh
bsub -J "sim_z[1-4]" < 02_simulate_z_LSF.sh
bsub -J "sim_magma[1-200]" < 03_magma_gene_LSF.sh
bsub -J "sim_merge[1-20]" < 03b_magma_merge_LSF.sh
bsub < 04_corpairs_once_LSF.sh
bsub -J "sim_catfish[1-20]" < 05_catfish_LSF.sh
bsub -J "sim_gsa[1-20]" < 06_magma_gsa_LSF.sh
bsub -J "sim_pxpval[1-20]" < 07a_pascalx_pval_LSF.sh       # R env
bsub -J "sim_pascalx[1-20]" < 07b_pascalx_LSF.sh           # PascalX env
bsub < 08_aggregate_LSF.sh
```
The pilot gives the run time of every step and the effect scale for the full
run (CATFISH omnibus power near 50% at the middle setting).

## Full run (MODE=full)
2 x 50 + 5 archetypes x 3 sizes x 50 = 850 replicates in 17 cells.
```
MODE=full EFFECT_SCALE=<from pilot> bsub < 01_make_scenarios_LSF.sh
MODE=full bsub -J "sim_z[1-17]" < 02_simulate_z_LSF.sh
MODE=full REP_FROM=1   bsub -J "sim_magma[1-1000]" < 03_magma_gene_LSF.sh   # replicates 1-100
MODE=full REP_FROM=101 bsub -J "sim_magma[1-1000]" < 03_magma_gene_LSF.sh   # 101-200 ... to 801
MODE=full bsub -J "sim_merge[1-850]" < 03b_magma_merge_LSF.sh
MODE=full bsub -J "sim_catfish[1-850]" < 05_catfish_LSF.sh
MODE=full bsub -J "sim_gsa[1-850]" < 06_magma_gsa_LSF.sh
MODE=full bsub -J "sim_pxpval[1-850]" < 07a_pascalx_pval_LSF.sh
MODE=full bsub -J "sim_pascalx[1-850]" < 07b_pascalx_LSF.sh
MODE=full bsub < 08_aggregate_LSF.sh
```
Step 04 is not rerun: the gene correlations depend only on the panel.

## Tested
The design, u, z, the z check and the aggregation were run end to end on a mock
panel (503 people, 3,000 SNPs) with a stand-in for plink2; 05_catfish.R was run
on mock MAGMA output. plink2's real --variant-score file layout is checked by
step 02 on the pilot (it stops if the columns are not as expected).
