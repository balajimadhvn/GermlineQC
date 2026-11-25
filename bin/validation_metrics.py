#!/usr/bin/env python3
"""
Turn per-sample vcfeval results into CAP/NABL-style validation tables:

  accuracy.tsv         TP / FP / FN, sensitivity, PPV, F1, per-base specificity (95% Wilson CI)
  concordance.tsv      pairwise genotype concordance for repeatability and reproducibility sets
  discordance.tsv      read depth / GQ of concordant vs discordant calls (root-cause analysis)
  titration.tsv        sensitivity / PPV by sequencing depth (limit of detection)
  acceptance.tsv       each acceptance criterion with observed value and PASS/FAIL
  samples.tsv          sample design + observed QC depth
"""
from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
from pathlib import Path

import yaml


# ---------- statistics ----------

def wilson(k: int, n: int, z: float = 1.96) -> tuple[float, float]:
    """95% Wilson score interval for a binomial proportion."""
    if n == 0:
        return (math.nan, math.nan)
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return (max(0.0, centre - half), min(1.0, centre + half))


def ratio(k: int, n: int) -> float:
    return k / n if n else math.nan


# ---------- VCF helpers ----------

def var_type(ref: str, alts: str) -> str:
    return "SNV" if all(len(a) == len(ref) == 1 for a in alts.split(",")) else "INDEL"


def count_by_type(vcf: Path) -> dict[str, int]:
    counts = {"SNV": 0, "INDEL": 0}
    if not vcf.exists():
        return counts
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] != "#":
                f = line.split("\t", 6)
                counts[var_type(f[3], f[4])] += 1
    return counts


def read_bed(path: Path) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = {}
    for line in path.read_text().splitlines():
        if line.strip():
            c, s, e = line.split("\t")[:3]
            regions.setdefault(c, []).append((int(s), int(e)))
    for c in regions:
        regions[c].sort()
    return regions


def bed_size(regions) -> int:
    return sum(e - s for iv in regions.values() for s, e in iv)


def in_bed(regions, chrom: str, pos0: int) -> bool:
    return any(s <= pos0 < e for s, e in regions.get(chrom, ()))


def pass_genotypes(vcfs: list[Path], regions) -> dict[tuple, tuple[str, int, int]]:
    """PASS, non-reference calls inside the evaluation region -> {(chrom,pos,ref,alt): (GT, DP, GQ)}."""
    calls = {}
    for vcf in vcfs:
        with gzip.open(vcf, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                if f[6] not in ("PASS", "."):
                    continue
                fmt = dict(zip(f[8].split(":"), f[9].split(":")))
                gt = fmt.get("GT", "./.").replace("|", "/")
                if gt in ("0/0", "./.", "."):
                    continue
                if in_bed(regions, f[0], int(f[1]) - 1):
                    dp = int(fmt["DP"]) if fmt.get("DP", ".").isdigit() else 0
                    gq = int(fmt["GQ"]) if fmt.get("GQ", ".").isdigit() else 0
                    calls[(f[0], int(f[1]), f[3], f[4])] = ("/".join(sorted(gt.split("/"))), dp, gq)
    return calls


# ---------- analyses ----------

def accuracy(sample: str, vcfeval_dir: Path) -> list[dict]:
    tpb = count_by_type(vcfeval_dir / "tp-baseline.vcf.gz")
    tpc = count_by_type(vcfeval_dir / "tp.vcf.gz")
    fn = count_by_type(vcfeval_dir / "fn.vcf.gz")
    fp = count_by_type(vcfeval_dir / "fp.vcf.gz")
    eval_bp = bed_size(read_bed(vcfeval_dir / "eval.bed"))
    rows = []
    for t in ("SNV", "INDEL"):
        sens = ratio(tpb[t], tpb[t] + fn[t])
        ppv = ratio(tpc[t], tpc[t] + fp[t])
        f1 = 2 * sens * ppv / (sens + ppv) if sens + ppv else math.nan
        s_lo, s_hi = wilson(tpb[t], tpb[t] + fn[t])
        p_lo, p_hi = wilson(tpc[t], tpc[t] + fp[t])
        rows.append(dict(sample=sample, type=t, tp=tpb[t], tp_call=tpc[t], fp=fp[t], fn=fn[t],
                         sensitivity=sens, sens_lo=s_lo, sens_hi=s_hi,
                         ppv=ppv, ppv_lo=p_lo, ppv_hi=p_hi, f1=f1))
    # per-base specificity: every confident base that is not a true variant is a negative
    fp_all = fp["SNV"] + fp["INDEL"]
    negatives = eval_bp - (tpb["SNV"] + tpb["INDEL"] + fn["SNV"] + fn["INDEL"])
    tn = negatives - fp_all
    lo, hi = wilson(tn, negatives)
    rows.append(dict(sample=sample, type="ALL", tp=tpb["SNV"] + tpb["INDEL"], fp=fp_all,
                     fn=fn["SNV"] + fn["INDEL"], tn=tn, eval_bp=eval_bp,
                     specificity=ratio(tn, negatives), spec_lo=lo, spec_hi=hi))
    return rows


def concordance(group: str, samples: list[str], vcfs: dict[str, list[Path]], regions):
    """Pairwise genotype concordance, plus depth/GQ of concordant vs discordant calls."""
    calls = {s: pass_genotypes(vcfs[s], regions) for s in samples}
    rows, quality = [], {"concordant": [], "discordant": []}
    for a, b in itertools.combinations(samples, 2):
        ca, cb = calls[a], calls[b]
        union = set(ca) | set(cb)
        both = set(ca) & set(cb)
        gt_match = sum(1 for k in both if ca[k][0] == cb[k][0])
        for k in union:
            x, y = ca.get(k), cb.get(k)
            if x and y and x[0] == y[0]:
                quality["concordant"].append(x)
            else:
                quality["discordant"].extend(v for v in (x, y) if v)
        lo, hi = wilson(gt_match, len(union))
        rows.append(dict(study=group, sample_a=a, sample_b=b, calls_a=len(ca), calls_b=len(cb),
                         shared=len(both), only_a=len(set(ca) - set(cb)), only_b=len(set(cb) - set(ca)),
                         genotype_match=gt_match, concordance=ratio(gt_match, len(union)),
                         conc_lo=lo, conc_hi=hi))
    summary = []
    for label, vals in quality.items():
        if vals:
            dps = sorted(v[1] for v in vals)
            summary.append(dict(study=group, calls=label, n=len(vals), median_dp=dps[len(dps) // 2],
                                frac_dp_lt15=sum(v[1] < 15 for v in vals) / len(vals),
                                frac_dp_lt20=sum(v[1] < 20 for v in vals) / len(vals),
                                frac_gq_lt20=sum(v[2] < 20 for v in vals) / len(vals)))
    return rows, summary


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", restval="")
        w.writeheader()
        for r in rows:
            w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v) for k, v in r.items()})


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcfeval", type=Path, nargs="+", required=True, help="vcfeval output directories")
    ap.add_argument("--qc", type=Path, nargs="+", required=True)
    ap.add_argument("--sheet", type=Path, required=True)
    ap.add_argument("--criteria", type=Path, required=True)
    ap.add_argument("--outdir", type=Path, required=True)
    a = ap.parse_args()
    a.outdir.mkdir(parents=True, exist_ok=True)

    design = {r["sample"]: r for r in csv.DictReader(a.sheet.open())}
    crit = yaml.safe_load(a.criteria.read_text())
    qc = {}
    for q in a.qc:
        d = json.loads(q.read_text())
        qc[d["sample"]] = {m["metric"]: m["value"] for m in d["metrics"]} | {"qc_status": d["overall"]}

    dirs = {d.name: d for d in a.vcfeval}
    acc = [row for s in sorted(dirs) for row in accuracy(s, dirs[s])]
    write_tsv(a.outdir / "accuracy.tsv", acc)

    # concordance uses the vcfeval output calls (already restricted to eval regions, PASS only)
    # vcfeval's call side is split into tp.vcf.gz + fp.vcf.gz (PASS calls inside the evaluation region)
    vcfs = {s: [dirs[s] / "tp.vcf.gz", dirs[s] / "fp.vcf.gz"] for s in dirs}
    conc, disc = [], []
    for group in ("repeatability", "reproducibility"):
        members = sorted(s for s in dirs if design.get(s, {}).get("study") == group)
        if len(members) > 1:
            regions = read_bed(dirs[members[0]] / "eval.bed")
            rows, summary = concordance(group, members, vcfs, regions)
            conc += rows
            disc += summary
    write_tsv(a.outdir / "concordance.tsv", conc)
    write_tsv(a.outdir / "discordance.tsv", disc)

    titr = []
    for r in acc:
        s = r["sample"]
        if design.get(s, {}).get("study") == "titration" and r["type"] != "ALL":
            titr.append(dict(sample=s, nominal_depth=int(design[s]["level"]),
                             observed_depth=qc.get(s, {}).get("mean_target_depth"),
                             type=r["type"], sensitivity=r["sensitivity"], ppv=r["ppv"], f1=r["f1"]))
    titr.sort(key=lambda x: (x["type"], x["nominal_depth"]))
    write_tsv(a.outdir / "titration.tsv", titr)

    samples = []
    for s, d in design.items():
        samples.append(dict(sample=s, study=d["study"], reference_material=d["reference_material"],
                            source=d["source"], nominal_depth=d["nominal_depth"],
                            observed_depth=qc.get(s, {}).get("mean_target_depth"),
                            pct_target_20x=qc.get(s, {}).get("pct_target_20x"),
                            qc_status=qc.get(s, {}).get("qc_status")))
    write_tsv(a.outdir / "samples.tsv", samples)

    # ---- acceptance criteria ----
    acc_set = [r for r in acc if design.get(r["sample"], {}).get("study") in ("accuracy", "repeatability", "reproducibility")]
    verdict = []

    def add(criterion, observed, threshold, detail):
        verdict.append(dict(criterion=criterion, observed=observed, threshold=threshold,
                            result="PASS" if observed >= threshold else "FAIL", detail=detail))

    for t in ("SNV", "INDEL"):
        rows = [r for r in acc_set if r["type"] == t]
        tp = sum(r["tp"] for r in rows); fn = sum(r["fn"] for r in rows); fp = sum(r["fp"] for r in rows)
        sens = ratio(tp, tp + fn)
        lo, _ = wilson(tp, tp + fn)
        add(f"{t} sensitivity (pooled)", sens, crit["accuracy"][t]["sensitivity"],
            f"{tp} TP / {tp + fn} truth variants across {len(rows)} runs; 95% CI lower bound {lo:.4f}")
        tpc = sum(r["tp_call"] for r in rows)
        ppv = ratio(tpc, tpc + fp)
        lo, _ = wilson(tpc, tpc + fp)
        add(f"{t} PPV (pooled)", ppv, crit["accuracy"][t]["ppv"],
            f"{fp} FP among {tpc + fp} PASS calls; 95% CI lower bound {lo:.4f}")

    spec_rows = [r for r in acc_set if r["type"] == "ALL"]
    tn = sum(r["tn"] for r in spec_rows); fp = sum(r["fp"] for r in spec_rows)
    add("Specificity (per base)", ratio(tn, tn + fp), crit["specificity"], f"{fp} FP over {tn + fp:,} reference bases")

    for group in ("repeatability", "reproducibility"):
        rows = [r for r in conc if r["study"] == group]
        if rows:
            worst = min(rows, key=lambda r: r["concordance"])
            add(f"{group.capitalize()} (min pairwise concordance)", worst["concordance"], crit[group],
                f"{len(rows)} pairs; worst {worst['sample_a']} vs {worst['sample_b']}")

    snv_titr = [r for r in titr if r["type"] == "SNV"]
    passing = [r["nominal_depth"] for r in snv_titr if r["sensitivity"] >= crit["lod_sensitivity"]]
    # LoD = lowest depth from which every higher depth also passes
    lod = None
    for r in snv_titr:
        if all(x["sensitivity"] >= crit["lod_sensitivity"] for x in snv_titr if x["nominal_depth"] >= r["nominal_depth"]):
            lod = r["nominal_depth"]
            break
    verdict.append(dict(criterion="Minimum depth for SNV sensitivity >= %.2f" % crit["lod_sensitivity"],
                        observed=lod if lod is not None else "not reached", threshold="report",
                        result="INFO", detail=f"depths meeting target: {passing}"))
    write_tsv(a.outdir / "acceptance.tsv", verdict)

    for v in verdict:
        obs = f"{v['observed']:.5f}" if isinstance(v["observed"], float) else v["observed"]
        print(f"{v['result']:<5} {v['criterion']:<48} {obs}")


if __name__ == "__main__":
    main()
