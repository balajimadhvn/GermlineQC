#!/usr/bin/env python3
"""
Collect per-sample QC metrics from fastp, samtools, Picard/GATK, mosdepth and bcftools,
then check each one against the lab's acceptance criteria (assets/qc_thresholds.yaml).
"""
from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import yaml


def fastp_metrics(path: Path) -> dict:
    d = json.loads(path.read_text())
    before, after = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    return {
        "raw_read_pairs": before["total_reads"] // 2,
        "q30_rate": round(after["q30_rate"], 4),
        "reads_passing_filter": round(after["total_reads"] / before["total_reads"], 4),
        "insert_size_peak": d.get("insert_size", {}).get("peak"),
    }


def flagstat_metrics(path: Path) -> dict:
    d = json.loads(path.read_text())["QC-passed reads"]
    return {"mapped_rate": round(d["mapped %"] / 100, 4),
            "properly_paired_rate": round(d["properly paired %"] / 100, 4)}


def dup_metrics(path: Path) -> dict:
    lines = path.read_text().splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith("LIBRARY"))
    row = dict(zip(lines[i].split("\t"), lines[i + 1].split("\t")))
    return {"duplication_rate": round(float(row["PERCENT_DUPLICATION"]), 4)}


def coverage_metrics(summary: Path, thresholds: Path) -> dict:
    mean = None
    for line in summary.read_text().splitlines()[1:]:
        f = line.split("\t")
        if f[0] == "total_region":
            mean = float(f[3])
    total, cov = 0, {}
    with gzip.open(thresholds, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        cols = header[4:]
        for line in fh:
            f = line.rstrip("\n").split("\t")
            size = int(f[2]) - int(f[1])
            total += size
            for c, v in zip(cols, f[4:]):
                cov[c] = cov.get(c, 0) + int(v)
    out = {"mean_target_depth": round(mean, 1) if mean is not None else None}
    for c, v in cov.items():
        out[f"pct_target_{c.lower()}"] = round(v / total, 4)
    return out


def variant_metrics(path: Path) -> dict:
    out = {}
    for line in path.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "SN":
            key = f[2].rstrip(":")
            if key == "number of SNPs":
                out["snvs"] = int(f[3])
            elif key == "number of indels":
                out["indels"] = int(f[3])
        elif f[0] == "TSTV" and len(f) > 4:
            out["ts_tv"] = float(f[4])
    return out


def evaluate(metrics: dict, criteria: dict) -> list[dict]:
    rows = []
    for name, value in metrics.items():
        rule = criteria.get(name, {})
        status = "INFO"
        if rule and value is not None:
            ok = True
            if "min" in rule:
                ok &= value >= rule["min"]
            if "max" in rule:
                ok &= value <= rule["max"]
            status = "PASS" if ok else ("WARN" if rule.get("level") == "warn" else "FAIL")
        rows.append({"metric": name, "value": value, "status": status,
                     "min": rule.get("min"), "max": rule.get("max")})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    for opt in ("fastp", "flagstat", "dup", "mosdepth", "thresholds", "bcfstats", "criteria", "json", "tsv"):
        ap.add_argument(f"--{opt}", type=Path, required=True)
    a = ap.parse_args()

    metrics = {}
    metrics |= fastp_metrics(a.fastp)
    metrics |= flagstat_metrics(a.flagstat)
    metrics |= dup_metrics(a.dup)
    metrics |= coverage_metrics(a.mosdepth, a.thresholds)
    metrics |= variant_metrics(a.bcfstats)

    rows = evaluate(metrics, yaml.safe_load(a.criteria.read_text()))
    overall = "FAIL" if any(r["status"] == "FAIL" for r in rows) else "PASS"
    a.json.write_text(json.dumps({"sample": a.sample, "overall": overall, "metrics": rows}, indent=2))
    with a.tsv.open("w") as fh:
        fh.write("sample\tmetric\tvalue\tmin\tmax\tstatus\n")
        for r in rows:
            fh.write(f"{a.sample}\t{r['metric']}\t{r['value']}\t{r['min'] or ''}\t{r['max'] or ''}\t{r['status']}\n")
    print(f"{a.sample}: {overall}")


if __name__ == "__main__":
    main()
