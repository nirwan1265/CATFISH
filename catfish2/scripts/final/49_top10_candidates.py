#!/usr/bin/env python3
"""
49_top10_candidates.py -- the ten highest-ranked CATFISH candidate genes of each
of the nine benchmark traits (gene score: GWAS + MAGMA + pathway layers), with
each gene's rank in every layer, its best pathway, and whether it is an Open
Targets gold-standard gene of the trait. Genes in the extended MHC (chr6
25-34 Mb, GRCh37) are skipped, as in the benchmark (45_benchmark_topk_all.R).
For reviewing which top candidates are known and which may be new.

REFS_FILE (default: 49_candidate_refs.csv next to this script): literature
evidence per trait and gene, curated by hand (columns trait, gene, evidence,
note, refs). evidence is "gene" (published genetic or functional evidence points
to this gene), "locus" (the gene lies at a locus reported for the trait; the
causal gene is unresolved) or "neighbour" (evidence points to a neighbouring
gene); refs are refs.bib keys separated by ";". Genes missing from the file are
reported and left blank, so rerun after editing it if the top lists change.

GWAS_CATALOG (optional): path to the GWAS Catalog associations file with
ontology annotations (.tsv, or the .zip it comes in; "All associations v1.0.2 -
with added ontology annotations" from https://www.ebi.ac.uk/gwas/docs/file-downloads).
If set, each gene gets the GWAS Catalog studies of the same trait (matched by
EFO/MONDO/HP term in MAPPED_TRAIT_URI, the terms of 45_benchmark_topk_all.R) that
list it as a reported or mapped gene: their number, and first author, year and
PubMed ID of the three earliest.
Reads  catfish2/results/06_candidate_genes_{ldl,hdl}.csv,
       catfish2/results/benchmark/candidates/06_candidate_genes_<trait>.csv,
       catfish2/results/45_gold_genes.csv, <CATFISH_DATA>/NCBI37.3.gene.loc (MAGMA gene locations)
Writes catfish2/results/49_top10_candidates.csv, catfish2/results/49_unresolved_loci.csv,
       catfish2/tex/auto/tab_top10_candidates.tex (Supplementary Table)
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

TERMS = {"ldl": ["EFO_0004611", "HP_0003124"], "hdl": ["EFO_0004612"], "cad": ["EFO_0001645", "EFO_0000378"],
         "brca": ["EFO_0000305"], "af": ["EFO_0000275"], "cd": ["EFO_0000384"], "uc": ["EFO_0000729"],
         "sle": ["EFO_0002690"], "asthma": ["MONDO_0004979", "EFO_0000270"]}
CAT = os.environ.get("GWAS_CATALOG", "")
cat = None
if CAT:
    import re
    cat = pd.read_csv(CAT, sep="\t", dtype=str, low_memory=False)   # pandas reads a single-file .zip too
    cat.columns = [c.strip() for c in cat.columns]
    need = ["PUBMEDID", "FIRST AUTHOR", "DATE", "REPORTED GENE(S)", "MAPPED_GENE", "MAPPED_TRAIT_URI"]
    miss = [c for c in need if c not in cat.columns]
    if miss: raise SystemExit("GWAS Catalog file lacks columns: %s" % miss)
    cat = cat[need].fillna("")
    split = lambda x: {g.strip() for g in re.split(r"[,;]| - | x ", x) if g.strip() and g.strip() not in ("NR", "intergenic")}
    cat["genes"] = [split(a) | split(b) for a, b in zip(cat["REPORTED GENE(S)"], cat["MAPPED_GENE"])]
    print("GWAS Catalog:", len(cat), "associations from", CAT)

def catalog_refs(t, gene):
    if cat is None: return {}
    terms = TERMS[t]
    m = cat[cat["MAPPED_TRAIT_URI"].apply(lambda u: any(x in u for x in terms)) & cat["genes"].apply(lambda s: gene in s)]
    st = m.drop_duplicates("PUBMEDID").sort_values("DATE")
    refs = ["%s %s (PMID %s)" % (r["FIRST AUTHOR"], r["DATE"][:4], r["PUBMEDID"]) for _, r in st.head(3).iterrows()]
    return dict(n_catalog_studies=len(st), catalog_refs="; ".join(refs))

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
    g = g[~g["mhc"]].copy()
    g["rank_gwas"] = g["rank_gwas"].rank(method="min")     # ranks among non-MHC genes, as in 50_lifted_demoted.py
    g["rank_magma"] = g["rank_magma"].rank(method="min")
    g = g.head(TOPN)
    for i, r in enumerate(g.itertuples(), 1):
        rows.append(dict(trait=lab, top=i, gene=r.symbol, rank_catfish=int(r.rank),
                         rank_gwas=r.rank_gwas, rank_magma=r.rank_magma, rank_pathway=r.rank_pathway,
                         layers=int(r.layers), best_pathway=pretty(r.best_pathway),
                         gold_standard=r.symbol in gold.get(t, set()),
                         chr=chr_.loc[r.Index], start=st.loc[r.Index], **catalog_refs(t, r.symbol)))
out = pd.DataFrame(rows)

REFS = os.environ.get("REFS_FILE", os.path.join(os.path.dirname(os.path.abspath(__file__)), "49_candidate_refs.csv"))
key = {v: k for k, v in TRAITS.items()}
if os.path.exists(REFS):
    ref = pd.read_csv(REFS, dtype=str).fillna("")
    ref["trait"] = ref["trait"].map(TRAITS)
    out = out.merge(ref.rename(columns={"note": "evidence_note"}), on=["trait", "gene"], how="left")
    miss = out.loc[out["evidence"].isna(), ["trait", "gene"]]
    if len(miss): print("No literature entry in", REFS, "for:\n", miss.to_string(index=False))
    out[["evidence", "evidence_note", "refs"]] = out[["evidence", "evidence_note", "refs"]].fillna("")
    print("evidence:", out["evidence"].value_counts().to_dict())
else:
    print("REFS_FILE not found:", REFS, "- table written without literature evidence")
    out["evidence"] = out["evidence_note"] = out["refs"] = ""
fo = os.path.join(RES, "49_top10_candidates.csv")
out.to_csv(fo, index=False)
## Genes at unresolved loci: the candidates for colocalization with expression QTL
fu = os.path.join(RES, "49_unresolved_loci.csv")
out[out["evidence"] == "locus"][["trait", "top", "gene", "chr", "start", "layers", "evidence_note", "refs"]].to_csv(fu, index=False)

## Supplementary table: one row per gene, split into floats of two traits.
## Gold-standard genes in bold; dagger = no pathway evidence (GWAS and MAGMA layers only).
import re
TEX = os.environ.get("TEX_AUTO", os.path.join(BASE, "catfish2/tex/auto"))
STOP = {"LDL", "HDL", "GWAS", "CAD", "AF", "IBD", "IBD5", "SLE", "UC", "CRISPR", "NADPH", "NF", "SK3", "IL", "BI", "DNA"}
def italic_genes(x):
    x = x.replace("_", "\\_").replace("NF-kB", "NF-$\\kappa$B")
    def f(m):
        w = m.group(0)
        if w in STOP or re.fullmatch(r"[A-Z]\d+[A-Z]", w): return w
        return "\\textit{%s}" % w
    x = re.sub(r"(?<![\w$\\])[A-Z][A-Z0-9]+(?![\w(])", f, x)
    return x.replace("}-\\textit{", "}--\\textit{").replace("}/\\textit{", "}/\\textit{")
def glab(r):
    s = "\\textit{%s}" % r["gene"].replace("_", "\\_")
    if r["gold_standard"]: s = "\\textbf{%s}" % s
    if r["layers"] < 3: s += "$^\\dagger$"
    return s
fmt = lambda x: "--" if pd.isna(x) else ("%d" % x if float(x).is_integer() else "%.1f" % x)
CLASS = {"gene": "Causal gene", "locus": "Unresolved locus", "neighbour": "Neighbour", "": ""}
cap = ("The ten highest-ranked CATFISH candidate genes of each of the nine human traits (gene score with the GWAS, "
       "MAGMA and pathway layers), with the published evidence linking each to the trait. "
       "Bold, Open Targets gold-standard gene of the trait; $^\\dagger$, no pathway evidence (ranked on the GWAS and MAGMA "
       "layers only). Evidence: Causal gene, published genetic or functional evidence points to this gene; Unresolved locus, the gene "
       "lies at a locus reported for the trait but the causal gene there has not been identified; Neighbour, published evidence points to a "
       "neighbouring gene at the same locus. GWAS and MAGMA rank, the gene's rank among all genes on each layer alone. Pathway ranks and best pathways are in \\texttt{49\\_top10\\_candidates.csv}.")
traits = list(dict.fromkeys(out["trait"]))
lines = []
for k, chunk in enumerate([traits[i:i + 2] for i in range(0, len(traits), 2)]):
    lines += ["\\begin{table*}[!p]", "\\centering",
              "" if k == 0 else "\\addtocounter{table}{-1}",   # continuation keeps the table number
              "\\caption{%s}" % (cap if k == 0 else "Top-ranked candidate genes and published evidence (continued)."),
              "\\label{tab:top10}" if k == 0 else "", "\\scriptsize",
              "\\begin{tabular}{p{0.11\\textwidth}p{0.03\\textwidth}p{0.09\\textwidth}r@{\\hspace{6pt}}r@{\\hspace{6pt}}p{0.08\\textwidth}>{\\raggedright\\arraybackslash}p{0.27\\textwidth}>{\\raggedright\\arraybackslash}p{0.16\\textwidth}}",
              "\\hline", "\\textbf{Trait} & \\textbf{Rank} & \\textbf{Gene} & \\textbf{GWAS rank} & \\textbf{MAGMA rank} & \\textbf{Evidence} & \\textbf{Published finding} & \\textbf{References} \\\\", "\\hline"]
    for t in chunk:
        g = out[out["trait"] == t]
        for j, (_, r) in enumerate(g.iterrows()):
            cite = "\\citealp{%s}" % ",".join(x.strip() for x in r["refs"].split(";") if x.strip()) if r["refs"] else ""
            lines.append("%s & %d & %s & %s & %s & %s & %s & %s \\\\" % (t if j == 0 else "", r["top"], glab(r), fmt(r["rank_gwas"]), fmt(r["rank_magma"]),
                         CLASS.get(r["evidence"], r["evidence"]), italic_genes(r["evidence_note"][:1].upper() + r["evidence_note"][1:]), cite))
        lines.append("\\hline")
    lines += ["\\end{tabular}", "\\end{table*}", ""]
ft = os.path.join(TEX, "tab_top10_candidates.tex")
open(ft, "w").write("\n".join(l for l in lines if l != "") + "\n")
print("wrote", fo, "and", ft)
