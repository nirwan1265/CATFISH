#!/usr/bin/env python3
"""
40_gold_standard_lipids.py -- LDL-C and HDL-C locus-to-gene gold standard for
the gene-prioritisation benchmark, from the Open Targets Genetics gold
standards (opentargets-archive/genetics-gold-standards, gold_standards/processed):
  gwas_gold_standards.191108.tsv   2019 set, GRCh37 and GRCh38 positions
  otg_gs_230511.json               2023 set (used for L2G), GRCh38 only; lifted
                                   to GRCh37 with the UCSC hg38ToHg19 chain
Traits: EFO_0004611 (LDL cholesterol measurement), EFO_0004612 (HDL cholesterol
measurement). Hypercholesterolaemia (HP_0003124) and total cholesterol are left
out: their loci list several genes per variant or mix the two traits.
One row per trait x gene (a gene listed by several studies is kept once, at
its highest confidence, with 2019 GRCh37 coordinates where they exist).

Reads  catfish2/data/gold_standards/{gwas_gold_standards.191108.tsv,
       otg_gs_230511.json, hg38ToHg19.over.chain.gz}
Writes catfish2/data/gold_standards/lipid_gold_standard.tsv
Usage (from catfish2/scripts/final):  python3 40_gold_standard_lipids.py
Needs: pip install pyliftover
"""
import csv, json, os, sys
from pyliftover import LiftOver

BASE = os.environ.get("CATFISH_BASE", os.path.expanduser("~/Documents/Github/CATFISH"))
GS = os.environ.get("GS_DIR", os.path.join(BASE, "catfish2/data/gold_standards"))
TRAITS = {"EFO_0004611": "ldl", "EFO_0004612": "hdl"}
RANK = {"High": 3, "Medium": 2, "Low": 1}

rows = []
with open(os.path.join(GS, "gwas_gold_standards.191108.tsv")) as fh:
    for x in csv.DictReader(fh, delimiter="\t"):
        ont = x["trait_info.ontology"].split("|")
        tr = [TRAITS[o] for o in ont if o in TRAITS]
        if not tr or not x["sentinel_variant.locus_GRCh37.position"]:
            continue
        rows.append(dict(trait=tr[0], gene=x["gold_standard_info.gene_id"],
                         chr=x["sentinel_variant.locus_GRCh37.chromosome"],
                         pos37=int(float(x["sentinel_variant.locus_GRCh37.position"])),
                         rsid=x["sentinel_variant.rsid"],
                         confidence=x["gold_standard_info.highest_confidence"],
                         study=x["association_info.otg_id"], source="2019"))

chain = os.path.join(GS, "hg38ToHg19.over.chain.gz")
n_fail = 0
if not os.path.exists(chain):
    sys.exit("missing " + chain + " (UCSC hg38ToHg19 liftover chain)")
lo = LiftOver(chain)
with open(os.path.join(GS, "otg_gs_230511.json")) as fh:
    for line in fh:
        if not line.strip():
            continue
        j = json.loads(line)
        tr = [TRAITS[o] for o in j["trait_info"]["ontology"] if o in TRAITS]
        if not tr:
            continue
        loc = j["sentinel_variant"]["locus_GRCh38"]
        hit = lo.convert_coordinate("chr" + str(loc["chromosome"]), int(loc["position"]) - 1)
        if not hit:
            n_fail += 1
            continue
        c37, p37 = hit[0][0].replace("chr", ""), hit[0][1] + 1
        rows.append(dict(trait=tr[0], gene=j["gold_standard_info"]["gene_id"], chr=c37, pos37=p37,
                         rsid="", confidence=j["gold_standard_info"]["highest_confidence"],
                         study=j["association_info"]["otg_id"], source="2023"))

## one row per trait x gene: highest confidence, then 2019 coordinates first
best = {}
for r in rows:
    k = (r["trait"], r["gene"])
    key = (RANK.get(r["confidence"], 0), r["source"] == "2019")
    if k not in best or key > best[k][0]:
        best[k] = (key, r)
out = sorted((v[1] for v in best.values()), key=lambda r: (r["trait"], int(r["chr"]) if r["chr"].isdigit() else 99, r["pos37"]))
cols = ["trait", "gene", "chr", "pos37", "rsid", "confidence", "study", "source"]
with open(os.path.join(GS, "lipid_gold_standard.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in out:
        fh.write("\t".join(str(r[c]) for c in cols) + "\n")
n = {t: sum(1 for r in out if r["trait"] == t) for t in ("ldl", "hdl")}
print(f"input rows {len(rows)} (liftover failures {n_fail}); kept LDL {n['ldl']}, HDL {n['hdl']} trait-gene pairs")
