#!/usr/bin/env python3

from pathlib import Path
import math
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import beta, binomtest, chi2


def lambda_of(pvals):
    p = pd.Series(pvals)
    p = p[p.notna() & (p > 0) & (p <= 1)]
    if len(p) == 0:
      return float("nan")
    chisq = chi2.isf(p.to_numpy(), 1)
    return float(np.median(chisq) / chi2.ppf(0.5, 1))


def read_inputs(in_dir):
    files = sorted(Path(in_dir).glob("pathway_pvals_perm_*.csv"))
    if not files:
        raise SystemExit(f"no per-perm CSVs in {in_dir}")
    dat = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    return files, dat


def build_type1(dat, pcols):
    rows = []
    for col in pcols:
        p = dat[col]
        p = p[p.notna() & np.isfinite(p) & (p >= 0) & (p <= 1)]
        n = len(p)
        if n < 1:
            continue
        for alpha in (0.05, 0.01):
            k = int((p <= alpha).sum())
            bt = binomtest(k, n, alpha)
            ci_lo, ci_hi = bt.proportion_ci(confidence_level=0.95)
            rows.append({
                "method": col,
                "alpha": alpha,
                "n": n,
                "observed": k / n,
                "ci_lo": ci_lo,
                "ci_hi": ci_hi,
                "calibrated": ci_lo <= alpha <= ci_hi,
            })
    return pd.DataFrame(rows)


def build_lambda(dat, pcols):
    return pd.DataFrame({
        "method": pcols,
        "lambda": [lambda_of(dat[col]) for col in pcols],
    }).dropna()


def qq_dataframe(pvals):
    p = pd.Series(pvals)
    p = p[p.notna() & np.isfinite(p) & (p > 0) & (p <= 1)].sort_values().to_numpy()
    n = len(p)
    probs = (np.arange(1, n + 1) - 0.5) / n
    return pd.DataFrame({
        "expected": -np.log10(probs),
        "observed": -np.log10(p),
    })


def qq_band(n):
    k = np.arange(1, n + 1)
    probs = (k - 0.5) / n
    return pd.DataFrame({
        "expected": -np.log10(probs),
        "lo": -np.log10(beta.ppf(0.975, k, n - k + 1)),
        "hi": -np.log10(beta.ppf(0.025, k, n - k + 1)),
    })


def save_plots(dat, head_col, n_files, out_dir):
    qd = qq_dataframe(dat[head_col])
    band = qq_band(len(qd))
    lam = lambda_of(dat[head_col])

    plt.figure(figsize=(6, 6), dpi=200)
    plt.fill_between(band["expected"], band["lo"], band["hi"], color="lightgray", alpha=0.6)
    plt.plot([qd["expected"].min(), qd["expected"].max()], [qd["expected"].min(), qd["expected"].max()],
             linestyle=(0, (4, 4)), color="red")
    plt.scatter(qd["expected"], qd["observed"], s=6, alpha=0.5, c="black")
    plt.xlabel("Expected − log10(p)")
    plt.ylabel("Observed − log10(p)")
    plt.title("CATFISH omnibus under phenotype-permutation null", fontweight="bold")
    plt.suptitle(f"{n_files} permutations, {len(qd)} pooled pathway p-values; lambda = {lam:.3f}",
                 y=0.92, fontsize=10)
    plt.tight_layout()
    plt.savefig(Path(out_dir) / "QQ_omnibus_null.png")
    plt.close()

    p = pd.Series(dat[head_col])
    p = p[p.notna() & np.isfinite(p) & (p >= 0) & (p <= 1)]
    bins = np.linspace(0, 1, 21)
    expected = len(p) / 20

    plt.figure(figsize=(6, 4), dpi=200)
    plt.hist(p, bins=bins, color="steelblue", edgecolor="white")
    plt.axhline(expected, color="red", linestyle=(0, (4, 4)))
    plt.xlabel("Omnibus p-value (null)")
    plt.ylabel("Count")
    plt.title("Null p-values should be ~Uniform(0,1)")
    plt.tight_layout()
    plt.savefig(Path(out_dir) / "hist_omnibus_null.png")
    plt.close()


def main():
    if len(sys.argv) < 3:
        raise SystemExit("Usage: 05_aggregate_and_plot_py.py <dir_with_perm_csvs> <output_dir>")

    in_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    files, dat = read_inputs(in_dir)
    pcols = [col for col in [
        "omni_p_final", "omni_p_mvn", "omni_p_analytic",
        "acat_p", "fisher_p", "minp_p_analytic",
        "stouffer_p_analytic", "tfisher_p_analytic",
        "acat_p_mvn_cal", "fisher_p_mvn_cal", "tfisher_p_mvn_cal",
        "minp_p_mvn_cal", "stouffer_p_mvn_cal", "omni_p_mvn_compcal",
    ] if col in dat.columns]

    type1 = build_type1(dat, pcols)
    lam = build_lambda(dat, pcols)
    type1.to_csv(out_dir / "type1_error_table.csv", index=False)
    lam.to_csv(out_dir / "lambda_table.csv", index=False)

    head_col = "omni_p_final" if "omni_p_final" in pcols else pcols[0]
    save_plots(dat, head_col, len(files), out_dir)

    print(f"[agg-py] pooling {len(files)} permutation files")
    print("\n=== Empirical type-I error (target = alpha) ===")
    print(type1.to_string(index=False))
    print("\n=== Inflation lambda (1.0 = perfectly calibrated) ===")
    print(lam.to_string(index=False))
    print(f"\n[agg-py] figures + tables written to {out_dir}")


if __name__ == "__main__":
    main()
