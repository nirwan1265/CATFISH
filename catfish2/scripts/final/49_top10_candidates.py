#!/usr/bin/env python3
"""
49_top10_candidates.py -- the ten highest-ranked CATFISH candidate genes of each
of the nine benchmark traits (gene score: GWAS + MAGMA + pathway layers), with
each gene's rank in every layer, its best pathway, and whether it is an Open
Targets gold-standard gene of the trait. Genes in the extended MHC (chr6
25-34 Mb, GRCh37) are skipped, as in the benchmark (45_benchmark_topk_all.R).
For reviewing which top candidates are known and which may be new.
Reads  catfish2/results/06_candidate_genes_{ldl,hdl}.csv,
       catfish2/results/benchmark/candidates/06_candidate_genes_<trait>.csv,
       catfish2/results/45_gold_genes.csv, <CATFISH_DATA>/NCBI37.3.gene.loc (MAGMA gene locations)
Writes catfish2/results/49_top10_candidates.csv
  python3 49_top10_candidates.py
"""
import os
import pandas as pd

BASE = os.environ.get("CATFISH_BASE", os.path.expanduser("~/Documents/Github/CATFISH"))
RES = os.path.join(BASE, "catfish2/results")
TOPN = int(os.environ.get("TOPN", "10"))
TRAITS = {"ldl": "LDL-C", "hdl": "HDL-C", "cad": "Coronary artery disease", "brca": "Breast cancer",
          "af": "Atrial fibrillation", "cd": "Crohn's disease", "uc": "Ulcerative colitis",
          "sle": "Systemic lupus erythematosus", "asthma": "Asthma"}

DATA = os.environ.get("CATFISH_DATA", os.path.expanduser("~/Documents/Github/MAGCAT_data"))
loc = pd.read_csv(os.path.join(DATA, "NCBI37.3.gene.loc"), sep=r"\s+", header=None, usecols=range(6),
                  names=["entrez", "CHR", "START", "END", "strand", "NAME"], dtype={"entrez": str, "CHR": str})
loc = loc.drop_duplicates("entrez").set_index("entrez")
gold = pd.read_csv(os.path.join(RES, "45_gold_genes.csv"))
gold = gold.groupby("trait")["symbol"].apply(set).to_dict()

def pretty(p):
    if not isinstance(p, str) or not p:
        return ""
    p = p.replace("REACTOME_", "").replace("_", " ").lower()
    return p[:1].upper() + p[1:]

rows = []
for t, lab in TRAITS.items():
    f = os.path.join(RES, f"06_candidate_genes_{t}.csv")
    if t not in ("ldl", "hdl") or not os.path.exists(f):
        f = os.path.join(RES, "benchmark", "candidates", f"06_candidate_genes_{t}.csv")
    g = pd.read_csv(f).sort_values("rank")
    g["gene"] = g["gene"].astype(str)
    chr_ = g["gene"].map(loc["CHR"]).astype(str)
    st = g["gene"].map(loc["START"]); en = g["gene"].map(loc["END"])
    g["mhc"] = (chr_ == "6") & (st <= 34e6) & (en >= 25e6)
    g = g[~g["mhc"]].head(TOPN)
    for i, r in enumerate(g.itertuples(), 1):
        rows.append(dict(trait=lab, top=i, gene=r.symbol, rank_catfish=int(r.rank),
                         rank_gwas=r.rank_gwas, rank_magma=r.rank_magma, rank_pathway=r.rank_pathway,
                         layers=int(r.layers), best_pathway=pretty(r.best_pathway),
                         gold_standard=r.symbol in gold.get(t, set()),
                         chr=chr_.loc[r.Index], start=st.loc[r.Index]))
out = pd.DataFrame(rows)
fo = os.path.join(RES, "49_top10_candidates.csv")
out.to_csv(fo, index=False)
print("wrote", fo)
