#!/usr/bin/env python3
"""
50_lifted_demoted.py -- which genes the CATFISH gene score moves into, and out
of, the genome-wide top N compared with ranking genes by their strongest GWAS
variant alone, for the nine benchmark traits, and how many of each are Open
Targets gold-standard genes of the trait.

  lifted : in the top N of the gene score but not of the GWAS layer
  dropped: in the top N of the GWAS layer but not of the gene score
A lifted gene is also marked "by pathway" when it is outside the top N of the
score without the pathway layer (S_GWAS + S_MAGMA), so the pathway layer, not
MAGMA, brought it in. Genes in the extended MHC (chr6 25-34 Mb, GRCh37) are
removed and all ranks recomputed among the remaining genes, as in
45_benchmark_topk_all.R and 49_top10_candidates.py.

TOPN (default 50): the top-N cutoff.
Reads  catfish2/results/06_candidate_genes_{ldl,hdl}.csv,
       catfish2/results/benchmark/candidates/06_candidate_genes_<trait>.csv,
       catfish2/results/45_gold_genes.csv, <CATFISH_DATA>/NCBI37.3.gene.loc
Writes catfish2/results/50_lifted_dropped_summary.csv,
       catfish2/results/50_lifted_dropped_genes.csv,
       catfish2/tex/auto/tab_lifted_dropped.tex       (Results table)
       catfish2/tex/auto/tab_lifted_genes.tex         (Supplementary table)
       catfish2/tex/auto/lifted.tex, lifted_stub.tex   (number macros for the text)
  python3 50_lifted_demoted.py
"""
import os
import pandas as pd

BASE = os.environ.get("CATFISH_BASE", os.path.expanduser("~/Documents/Github/CATFISH"))
DATA = os.environ.get("CATFISH_DATA", os.path.expanduser("~/Documents/Github/MAGCAT_data"))
RES = os.path.join(BASE, "catfish2/results")
TEX = os.environ.get("TEX_AUTO", os.path.join(BASE, "catfish2/tex/auto"))
N = int(os.environ.get("TOPN", "50"))
TRAITS = {"ldl": "LDL-C", "hdl": "HDL-C", "cad": "Coronary artery disease", "brca": "Breast cancer",
          "af": "Atrial fibrillation", "cd": "Crohn's disease", "uc": "Ulcerative colitis",
          "sle": "Systemic lupus erythematosus", "asthma": "Asthma"}

loc = pd.read_csv(os.path.join(DATA, "NCBI37.3.gene.loc"), sep=r"\s+", header=None, usecols=range(6),
                  names=["entrez", "CHR", "START", "END", "strand", "NAME"], dtype={"entrez": str, "CHR": str})
loc = loc.drop_duplicates("entrez").set_index("entrez")
gold = pd.read_csv(os.path.join(RES, "45_gold_genes.csv")).groupby("trait")["symbol"].apply(set).to_dict()

summ, genes = [], []
for t, lab in TRAITS.items():
    f = os.path.join(RES, f"06_candidate_genes_{t}.csv")
    if t not in ("ldl", "hdl") or not os.path.exists(f):
        f = os.path.join(RES, "benchmark", "candidates", f"06_candidate_genes_{t}.csv")
    g = pd.read_csv(f)
    g["gene"] = g["gene"].astype(str)
    c = g["gene"].map(loc["CHR"]).astype(str); s = g["gene"].map(loc["START"]); e = g["gene"].map(loc["END"])
    g = g[~((c == "6") & (s <= 34e6) & (e >= 25e6))].copy()
    lay = g[["S_GWAS", "S_MAGMA"]].fillna(0)
    g["r_cat"] = g["Score"].rank(ascending=False, method="first")
    g["r_gwas"] = g["rank_gwas"].rank(method="min")
    g["r_magma"] = g["rank_magma"].rank(method="min")
    g["r_nopath"] = (lay["S_GWAS"] + lay["S_MAGMA"]).rank(ascending=False, method="first")
    g["gold"] = g["symbol"].isin(gold.get(t, set()))
    top10 = g[g["r_cat"] <= 10]
    lifted = g[(g["r_cat"] <= N) & (g["r_gwas"] > N)].copy()
    lifted["by_pathway"] = lifted["r_nopath"] > N
    dropped = g[(g["r_gwas"] <= N) & (g["r_cat"] > N)]
    kept = g[(g["r_gwas"] <= N) & (g["r_cat"] <= N)]
    summ.append(dict(trait=lab, top10_gwas_min=int(top10["r_gwas"].min()), top10_gwas_max=int(top10["r_gwas"].max()),
                     top10_outside_gwas_top10=int((top10["r_gwas"] > 10).sum()),
                     lifted=len(lifted), lifted_by_pathway=int(lifted["by_pathway"].sum()),
                     lifted_gold=int(lifted["gold"].sum()), dropped=len(dropped), dropped_gold=int(dropped["gold"].sum()),
                     kept=len(kept), kept_gold=int(kept["gold"].sum()),
                     gold_total=len(gold.get(t, set()) & set(g["symbol"]))))
    for kind, d in (("lifted", lifted), ("dropped", dropped)):
        for r in d.sort_values("r_cat").itertuples():
            genes.append(dict(trait=lab, change=kind, gene=r.symbol, rank_gwas=int(r.r_gwas), rank_magma=int(r.r_magma),
                              rank_nopathway=int(r.r_nopath), rank_catfish=int(r.r_cat), gold_standard=r.gold,
                              by_pathway=bool(getattr(r, "by_pathway", False)), best_pathway=r.best_pathway if isinstance(r.best_pathway, str) else ""))
S = pd.DataFrame(summ); G = pd.DataFrame(genes)
S.to_csv(os.path.join(RES, "50_lifted_dropped_summary.csv"), index=False)
G.to_csv(os.path.join(RES, "50_lifted_dropped_genes.csv"), index=False)
tot = S[["top10_outside_gwas_top10", "lifted", "lifted_by_pathway", "lifted_gold", "dropped", "dropped_gold", "kept", "kept_gold"]].sum()
print(S.to_string(index=False)); print("totals:", tot.to_dict())

## Results table
L = ["\\begin{table*}[!t]", "\\centering",
     "\\caption{Genes the CATFISH gene score moves into and out of the genome-wide top %d, compared with ranking genes by "
     "their strongest GWAS variant, in nine human traits (extended MHC removed). Lifted, in the top %d of the score but not "
     "of the GWAS layer; by pathway, lifted genes also outside the top %d of the score without the pathway layer; dropped, "
     "in the top %d of the GWAS layer but not of the score; gold, Open Targets gold-standard genes of the trait. Top 10, range "
     "of GWAS-layer ranks of the ten highest-scoring genes. Lifted genes are listed in Supplementary Table~\\ref{tab:lifted_genes}.}" % (N, N, N, N),
     "\\label{tab:lifted_dropped}", "\\small", "\\begin{tabular}{lrrrrrr}", "\\hline",
     "& \\multicolumn{1}{c}{Top 10} & \\multicolumn{3}{c}{Lifted} & \\multicolumn{2}{c}{Dropped} \\\\",
     "Trait & GWAS ranks & All & By pathway & Gold & All & Gold \\\\", "\\hline"]
for r in S.itertuples():
    L.append("%s & %d--%d & %d & %d & %d & %d & %d \\\\" % (r.trait, r.top10_gwas_min, r.top10_gwas_max, r.lifted,
             r.lifted_by_pathway, r.lifted_gold, r.dropped, r.dropped_gold))
L += ["\\hline", "Total & & %d & %d & %d & %d & %d \\\\" % (tot.lifted, tot.lifted_by_pathway, tot.lifted_gold, tot.dropped, tot.dropped_gold),
      "\\hline", "\\end{tabular}", "\\end{table*}"]
open(os.path.join(TEX, "tab_lifted_dropped.tex"), "w").write("\n".join(L) + "\n")

## Supplementary table: lifted genes, GWAS rank -> score rank
def lab(r):
    s = "\\textit{%s}" % r["gene"].replace("_", "\\_")
    if r["gold_standard"]: s = "\\textbf{%s}" % s
    if r["by_pathway"]: s += "$^\\ast$"
    return "%s (%d$\\rightarrow$%d)" % (s, r["rank_gwas"], r["rank_catfish"])
L = ["\\begin{table*}[!t]", "\\centering",
     "\\caption{Genes lifted into the genome-wide top %d of the CATFISH gene score from outside the top %d of the GWAS layer, "
     "in score order, with each gene's GWAS-layer rank and score rank (extended MHC removed). Bold, Open Targets gold-standard "
     "gene of the trait; $^\\ast$, also outside the top %d of the score without the pathway layer. Ranks on every layer are in "
     "\\texttt{50\\_lifted\\_dropped\\_genes.csv}.}" % (N, N, N),
     "\\label{tab:lifted_genes}", "\\footnotesize", "\\begin{tabular}{p{0.16\\textwidth}>{\\raggedright\\arraybackslash}p{0.78\\textwidth}}", "\\hline",
     "\\textbf{Trait} & \\textbf{Lifted genes (GWAS rank $\\rightarrow$ score rank)} \\\\", "\\hline"]
for t in TRAITS.values():
    d = G[(G["trait"] == t) & (G["change"] == "lifted")]
    L.append("%s & %s \\\\" % (t, ", ".join(lab(r) for _, r in d.iterrows())))
L += ["\\hline", "\\end{tabular}", "\\end{table*}"]
open(os.path.join(TEX, "tab_lifted_genes.tex"), "w").write("\n".join(L) + "\n")

## Macros for the text
M = {"LiftN": N, "LiftTopTenOut": tot.top10_outside_gwas_top10, "LiftAll": tot.lifted, "LiftPath": tot.lifted_by_pathway,
     "LiftGold": tot.lifted_gold, "LiftPathGold": int((G[G["change"] == "lifted"]["gold_standard"] & G[G["change"] == "lifted"]["by_pathway"]).sum()), "DropAll": tot.dropped, "DropGold": tot.dropped_gold, "KeepGold": tot.kept_gold, "KeepAll": tot.kept}
M["LiftGoldPct"] = "%.1f" % (100 * tot.lifted_gold / tot.lifted); M["KeepGoldPct"] = "%.0f" % (100 * tot.kept_gold / tot.kept)
open(os.path.join(TEX, "lifted.tex"), "w").write("".join("\\renewcommand{\\n%s}{%s}\n" % (k, v) for k, v in M.items()))
open(os.path.join(TEX, "lifted_stub.tex"), "w").write("".join("\\providecommand{\\n%s}{\\textbf{??}}\n" % k for k in M))
print("wrote tables and macros to", TEX)
