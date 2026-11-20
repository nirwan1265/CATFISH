# catfish2/scripts/final

Scripts for the analyses added in Sections 3.10 and 3.11 of `tex/catfish_wy.tex` and the gene-score checks in Section 3.7. Run each one from this folder with `Rscript <script>`. Every input and output path comes from `CATFISH_BASE` (default `~/Documents/Github/CATFISH`), so the scripts do not depend on the folder they are run from.

| Script | Paper | Reads | Writes |
|---|---|---|---|
| `35_simulation_figures.R` | §3.10, Figure 18 and Table 17 (simulation, CATFISH vs MAGMA) | `results/simulation_check_catfish/run_full/summary.csv` | `figures/Fig22_simulation_catfish_magma.png`, `tex/auto/tab_simulation_compare.tex` |
| `36_prioritize_comp.R` | §3.11, gene score with a competitive or joint pathway layer | `results/gc_real_comp/` | `results/prioritize_comp/` (it runs `../06_candidate_genes.R` and `../25_sorghum_gene_prioritize.R`) |
| `37_prioritize_comp_recovery.R` | §3.11, established-gene recovery for each pathway layer | the outputs of 36 and the paper's 06/25 outputs | `results/prioritize_comp/37_recovery_*.csv` |
| `38_competitive_real.R` | §3.11, Table 18 (self-contained vs competitive in the real data) | `results/gc_real_comp/` | `tex/auto/tab_competitive_real.tex`, `results/38_competitive_real.csv` |
| `39_logo_recovery.R` | §2.5 and §3.7, Table `tab:logo_recovery` (layer and leave-one-gene-out checks) | the 06 outputs, `results/human_LDL/catfish/`, `results/gc_real_logo/` | `tex/auto/tab_logo_recovery.tex`, `results/39_logo_recovery_*.csv` |

## Cluster steps behind them (in `scripts/hpc/`)

These stay in `scripts/hpc/` because they run on Hazel.

- `simulation_check_catfish/`: the full simulation pipeline, steps 00–08 (see its `README.md`). It produces `summary.csv` for 35.
- `gc_real/wy_gc_real.R` with `COMP=1` (`wy_gc_real_LSF.sh`): CATFISH's competitive test on the real data. Its output `out_comp/` is copied to `results/gc_real_comp/` for 36 and 38.
- `gc_real/wy_logo.R` (`wy_logo_LSF.sh`): leave-one-gene-out pathway statistics. Its output `out_logo/` is copied to `results/gc_real_logo/` for 39.

## Order

1. **Simulation:** cluster steps 00–08, then `35`.
2. **Competitive test on the real data:** cluster `wy_gc_real` with `COMP=1`, rsync `out_comp/`, then `38`, then `36` and `37`.
3. **Leave one gene out:** cluster `wy_logo`, rsync `out_logo/`, then `39`.

`36` passes extra environment variables to `06` and `25` (`CANDIDATE_OUT_DIR`; `WY_DIR`, `OUT_DIR`, `PATHWAY_STAT_COL`). When those variables are not set, both scripts behave exactly as before. For sorghum, also set `PATHWAY_FILE` to the SorghumBicolorCyc `.SORBI` file, and for HDL-C EUR set `GWAS_FILE` to `MAGCAT_data/hdl_eur.neglog10.txt`.
4. **BAP stem composition and the chr4 sugar/fibre locus (Section 3.12, 4.8):** cluster `hpc/bap_gwas/01`–`05` (Brix_maturity, NDF_pct_DM, NFC_pct_DM, 3 PCs), rsync `catfish/<trait>/` to `results/bap_traits/<trait>/`, `prioritize/` and `chr4_locus/` to `results/bap_traits/`, then `43` (Figure 23) and `44` (numbers `tex/auto/bap_chr4.tex`, table `tex/auto/tab_bap_chr4.tex`).
