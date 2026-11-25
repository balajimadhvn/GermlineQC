#!/usr/bin/env python3
"""
Score acmg_classify.py output against ClinVar consensus classifications (>= 2 review stars).

ClinVar classes are collapsed to P/LP and B/LB; predictions are reported in five tiers and three
(P/LP, VUS, B/LB). The numbers that matter clinically:

  * concordance          prediction in the same direction as ClinVar
  * VUS rate             the classifier could not reach a decision (safe but unhelpful)
  * opposite-direction   ClinVar P/LP called B/LB or the reverse (the errors that harm patients)

Outputs: summary.tsv, confusion.tsv, by_consequence.tsv, criteria.tsv, discordant.tsv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import math
from collections import Counter
from pathlib import Path

TIER3 = {"Pathogenic": "P/LP", "Likely pathogenic": "P/LP", "Uncertain significance": "VUS",
         "Likely benign": "B/LB", "Benign": "B/LB"}
CLINVAR3 = {"Pathogenic": "P/LP", "Likely_pathogenic": "P/LP", "Pathogenic/Likely_pathogenic": "P/LP",
            "Benign": "B/LB", "Likely_benign": "B/LB", "Benign/Likely_benign": "B/LB"}
FIVE = ["Pathogenic", "Likely pathogenic", "Uncertain significance", "Likely benign", "Benign"]


def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, c - h), min(1.0, c + h))


def consequence_group(csq: str) -> str:
    c = set(csq.split("&"))
    if c & {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}:
        return "null (LoF)"
    if "missense_variant" in c:
        return "missense"
    if c & {"inframe_insertion", "inframe_deletion", "disruptive_inframe_insertion",
            "disruptive_inframe_deletion", "stop_lost"}:
        return "in-frame / stop-loss"
    if "synonymous_variant" in c:
        return "synonymous"
    if "splice_region_variant" in c:
        return "splice region"
    if c & {"intron_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
            "upstream_gene_variant", "downstream_gene_variant"}:
        return "intronic / UTR / flanking"
    return "other / unannotated"


def read_truth(vcf: Path) -> dict[tuple, str]:
    truth = {}
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.split("\t", 8)
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            if info.get("CLNSIG") in CLINVAR3:
                truth[(f[0], int(f[1]), f[3], f[4])] = CLINVAR3[info["CLNSIG"]]
    return truth


def write(path: Path, rows: list[dict]) -> None:
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) else v) for k, v in r.items()})


def rates(rows: list[dict], truth_class: str) -> dict:
    sub = [r for r in rows if r["truth"] == truth_class]
    n = len(sub)
    same = sum(r["pred3"] == truth_class for r in sub)
    vus = sum(r["pred3"] == "VUS" for r in sub)
    opp = n - same - vus
    lo, hi = wilson(same, n)
    olo, ohi = wilson(opp, n)
    return dict(truth=truth_class, n=n, concordant=same, vus=vus, opposite=opp,
                concordance=same / n if n else math.nan, conc_lo=lo, conc_hi=hi,
                vus_rate=vus / n if n else math.nan,
                opposite_rate=opp / n if n else math.nan, opp_lo=olo, opp_hi=ohi)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--classified", type=Path, required=True)
    ap.add_argument("--truth", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    truth = read_truth(a.truth)
    rows = []
    for r in csv.DictReader(a.classified.open(), delimiter="\t"):
        k = (r["chrom"], int(r["pos"]), r["ref"], r["alt"])
        if k in truth:
            rows.append(dict(r, truth=truth[k], pred3=TIER3[r["acmg_class"]],
                             group=consequence_group(r["consequence"])))

    summary = [rates(rows, "P/LP"), rates(rows, "B/LB")]
    write(a.outdir / "summary.tsv", summary)

    conf = Counter((r["truth"], r["acmg_class"]) for r in rows)
    write(a.outdir / "confusion.tsv",
          [dict(truth=t, predicted=p, n=conf[(t, p)]) for t in ("P/LP", "B/LB") for p in FIVE])

    by_csq = []
    for g in sorted({r["group"] for r in rows}):
        for t in ("P/LP", "B/LB"):
            sub = [r for r in rows if r["group"] == g]
            x = rates(sub, t)
            if x["n"]:
                by_csq.append(dict(consequence=g, **x))
    write(a.outdir / "by_consequence.tsv", by_csq)

    crit = Counter()
    for r in rows:
        for c in filter(None, r["criteria"].split(",")):
            crit[(c.split("_")[0], r["truth"])] += 1
    write(a.outdir / "criteria.tsv",
          [dict(criterion=c, truth=t, n=n) for (c, t), n in sorted(crit.items())])

    disc = [dict(chrom=r["chrom"], pos=r["pos"], ref=r["ref"], alt=r["alt"], gene=r["gene"],
                 consequence=r["consequence"], hgvs_p=r["hgvs_p"], gnomad_faf=r["gnomad_faf"], revel=r["revel"],
                 criteria=r["criteria"], predicted=r["acmg_class"], clinvar=r["truth"])
            for r in rows if r["pred3"] not in ("VUS", r["truth"])]
    write(a.outdir / "discordant.tsv", disc or [dict(note="no opposite-direction calls")])

    print(f"{len(rows)} benchmark variants scored")
    for s in summary:
        print(f"  ClinVar {s['truth']:<5} n={s['n']:<5} concordant {s['concordance']:.3f}  "
              f"VUS {s['vus_rate']:.3f}  opposite {s['opposite_rate']:.4f}")


if __name__ == "__main__":
    main()
