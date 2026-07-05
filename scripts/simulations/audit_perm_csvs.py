#!/usr/bin/env python3

import csv
import glob
import math
import os
import sys
from collections import Counter, defaultdict


def is_finite_number(value):
    if value is None:
        return False
    s = str(value).strip()
    if s == "" or s.lower() in {"na", "nan", "null"}:
        return False
    try:
        x = float(s)
    except ValueError:
        return False
    return math.isfinite(x)


def read_file_summary(path):
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames or []
        n_rows = 0
        finite_counts = Counter()
        for row in reader:
            n_rows += 1
            for col, val in row.items():
                if is_finite_number(val):
                    finite_counts[col] += 1
    return {
        "file": path,
        "basename": os.path.basename(path),
        "n_rows": n_rows,
        "columns": tuple(fieldnames),
        "finite_counts": finite_counts,
    }


def main():
    if len(sys.argv) != 2:
        print("Usage: audit_perm_csvs.py <dir_with_pathway_pvals_perm_csvs>", file=sys.stderr)
        sys.exit(1)

    in_dir = sys.argv[1]
    files = sorted(glob.glob(os.path.join(in_dir, "pathway_pvals_perm_*.csv")))
    if not files:
        print(f"No pathway_pvals_perm_*.csv files found in {in_dir}", file=sys.stderr)
        sys.exit(2)

    summaries = [read_file_summary(path) for path in files]

    print(f"[audit] directory: {in_dir}")
    print(f"[audit] files: {len(summaries)}")

    row_counts = Counter(s["n_rows"] for s in summaries)
    print(f"[audit] row-count distribution: {dict(sorted(row_counts.items()))}")

    column_layouts = Counter(s["columns"] for s in summaries)
    print(f"[audit] distinct column layouts: {len(column_layouts)}")
    for idx, (cols, n) in enumerate(column_layouts.items(), start=1):
        print(f"  layout {idx}: {n} files, {len(cols)} columns")
        print(f"    {list(cols)}")

    tracked_cols = [
        "omni_p_final",
        "omni_p_mvn",
        "omni_p_analytic",
        "omni_p_mvn_compcal",
        "acat_p_mvn_cal",
        "fisher_p_mvn_cal",
        "tfisher_p_mvn_cal",
        "minp_p_mvn_cal",
        "stouffer_p_mvn_cal",
    ]

    print("\n[audit] finite-value counts by file for calibration-sensitive columns")
    header = ["file", "n_rows"] + tracked_cols
    print("\t".join(header))
    for s in summaries:
        row = [s["basename"], str(s["n_rows"])]
        for col in tracked_cols:
            row.append(str(s["finite_counts"].get(col, 0)))
        print("\t".join(row))

    suspicious = []
    for col in tracked_cols:
        counts = Counter(s["finite_counts"].get(col, 0) for s in summaries)
        if len(counts) > 1:
            suspicious.append((col, counts))

    if suspicious:
        print("\n[audit] columns with inconsistent finite-count patterns across files")
        for col, counts in suspicious:
            print(f"  {col}: {dict(sorted(counts.items()))}")
    else:
        print("\n[audit] all tracked columns are consistent across files")

    by_pattern = defaultdict(list)
    for s in summaries:
        key = tuple((col, s["finite_counts"].get(col, 0)) for col in tracked_cols)
        by_pattern[key].append(s["basename"])

    if len(by_pattern) > 1:
        print("\n[audit] distinct file patterns detected")
        for idx, (pattern, members) in enumerate(by_pattern.items(), start=1):
            print(f"  pattern {idx}: {len(members)} files")
            print(f"    files: {', '.join(members)}")
            print(f"    counts: {dict(pattern)}")
    else:
        print("\n[audit] all files share the same tracked finite-count pattern")


if __name__ == "__main__":
    main()
