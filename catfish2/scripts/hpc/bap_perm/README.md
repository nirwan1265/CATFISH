# BAP phenotype-permutation null for the Westfall-Young paper

Dry tons per acre is shuffled among the BAP lines 100 times. Each shuffle goes
through the real pipeline: vcf2gwas LMM (3 PCs), MAGMA per chromosome, then
Westfall-Young with the settings of `23_westfall_young_sorghum.R`. With no true
association, the Westfall-Young family-wise error should stay at or below 5%.

Scripts live in `/share/maize/ntanduk/CATFISH/simulation/bap_perm_wy`; outputs
go to `run/` under it. Every job is `-q sara` with no memory request, and runs
one chromosome per job.

| Step | Submit | Jobs |
|---|---|---|
| 0 | `cd $CATFISH_REPO && git pull` (needs `catfish_wy_pathways()`) | - |
| 0 | `bsub < 00_prepare_genome_LSF.sh` (vcf2gwas env active): per-chromosome VCFs of the phenotyped lines; genome-wide kinship and 3 PCs, once | 1 |
| 1 | `bsub < 01_make_permutations_LSF.sh` | 1 |
| 2 | `bsub < 02_gwas_array_LSF.sh` (vcf2gwas env active): vcf2gwas per (permutation, chromosome) with the genome-wide `-k` kinship and `-cf` PCs | 1000 |
| 2b | `bsub < 02b_annotate_LSF.sh`: MAGMA annotation per chromosome, built once as the real analysis did (.bim SNP locations, sorghum.genes.loc, 25,25 kb) | 10 |
| 3 | `bsub < 03_magma_array_LSF.sh`: MAGMA per (permutation, chromosome) | 1000 |
| 4 | `bsub < 04_corpairs_once_LSF.sh` | 1 |
| 5 | `bsub < 05_wy_array_LSF.sh` | 100 |
| 6 | `bsub < 06_aggregate_LSF.sh` | 1 |

Submit each step only after the previous one has finished, and check it
before going on. The earlier per-chromosome GWAS version, which rebuilt the
kinship and PCs from one chromosome in every job, is kept in `_old_perchr/`;
the whole-genome version in `_old_wholegenome/`.

Settings, all in `config.sh`:
- **N_PERM** = 100; MASTER_SEED is the same as the earlier run.
- **GENE_MODEL** = `snp-wise=mean`, which gives the paper's gene p-value fast.
  `multi=snp-wise` also gives P_MULTI, at about 8 h per chromosome.
- **Annotation**: the real analysis's genotype-only per-chromosome `.genes.annot`
  (+/- 25 kb).
