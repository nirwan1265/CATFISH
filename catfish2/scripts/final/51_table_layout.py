#!/usr/bin/env python3
"""
51_table_layout.py -- page layout of the generated LaTeX tables, run after any
script that rewrites a table in tex/auto (35, 39, 43, 45, 49, 50, 07, ...).

  * every table spans the full width of its float: \\begin{tabular}{spec}
    becomes \\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}spec}
    (\\linewidth is the page width in a table* and the column width in a table)
  * supplementary tables may be placed where they are defined ([!htbp]), so in
    the one-column supplement each table follows its heading instead of being
    pushed to the top of the page; main-text tables keep their placement

The edit is idempotent: files already converted are left as they are.
Usage (from catfish2/scripts):  python3 final/51_table_layout.py
"""
import os, re

BASE = os.environ.get("CATFISH_BASE", os.path.expanduser("~/Documents/Github/CATFISH"))
AUTO = os.environ.get("TEX_AUTO", os.path.join(BASE, "catfish2/tex/auto"))
MAIN = {"tab_simulation_compare", "tab_benchmark_topk_all", "tab_lifted_dropped", "tab_bap_chr4"}
SUPP = {"tab_breakage", "tab_loo", "tab_components_cal", "tab_missing", "tab_lipid_pathways",
        "tab_magma_compare", "tab_competitive_real", "tab_lipid_candidates", "tab_replication_genes",
        "tab_apoe_locus", "tab_logo_recovery", "tab_logo_recovery_bench", "tab_pathway_gain_all",
        "tab_lifted_genes", "tab_top10_candidates", "tab_sorghum_null"}

def widen(s):
    s = re.sub(r"\\begin\{tabular\}\{", r"\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}", s)
    return s.replace("\\end{tabular}", "\\end{tabular*}")

def float_here(s):
    return re.sub(r"\\begin\{(table\*?)\}\[!?[tbp]+\]", r"\\begin{\1}[!htbp]", s)

for name in sorted(MAIN | SUPP):
    f = os.path.join(AUTO, name + ".tex")
    if not os.path.exists(f):
        print("missing", f); continue
    s0 = open(f).read(); s = widen(s0)
    if name in SUPP: s = float_here(s)
    if s != s0:
        open(f, "w").write(s); print("updated", name)
