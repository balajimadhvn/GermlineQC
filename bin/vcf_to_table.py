#!/usr/bin/env python3
"""Flatten a SnpEff + ClinVar annotated single-sample VCF into a tab-separated table."""
from __future__ import annotations

import gzip
import sys

IMPACT_RANK = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
COLUMNS = ["sample", "chrom", "pos", "ref", "alt", "var_type", "filter", "qual", "genotype",
           "depth", "gq", "vaf", "gene", "effect", "impact", "transcript", "hgvs_c", "hgvs_p",
           "clinvar_id", "clinvar_sig", "clinvar_review", "clinvar_disease"]


def parse_info(field: str) -> dict[str, str]:
    out = {}
    for item in field.split(";"):
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "true"
    return out


def top_annotation(ann: str) -> list[str]:
    """Pick the most severe SnpEff ANN entry (ties keep SnpEff's order)."""
    if not ann:
        return [""] * 16
    entries = [e.split("|") for e in ann.split(",")]
    entries.sort(key=lambda e: IMPACT_RANK.get(e[2], 9))
    return entries[0]


def genotype(gt: str) -> str:
    alleles = gt.replace("|", "/").split("/")
    if "." in alleles:
        return "nocall"
    return "hom_alt" if len(set(alleles)) == 1 and alleles[0] != "0" else (
        "hom_ref" if set(alleles) == {"0"} else "het")


def main(path: str, sample: str) -> None:
    w = sys.stdout.write
    w("\t".join(COLUMNS) + "\n")
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom, pos, _, ref, alt, qual, flt, info = f[:8]
            fmt = dict(zip(f[8].split(":"), f[9].split(":")))
            ad = [int(x) for x in fmt.get("AD", "").split(",") if x.isdigit()]
            vaf = f"{ad[1] / sum(ad):.3f}" if len(ad) > 1 and sum(ad) else ""
            i = parse_info(info)
            a = top_annotation(i.get("ANN", ""))
            var_type = "SNV" if len(ref) == len(alt) == 1 else "INDEL"
            row = [sample, chrom, pos, ref, alt, var_type, flt, qual, genotype(fmt.get("GT", "./.")),
                   fmt.get("DP", ""), fmt.get("GQ", ""), vaf,
                   a[3], a[1], a[2], a[6], a[9], a[10],
                   i.get("CLNID", ""), i.get("CLNSIG", ""), i.get("CLNREVSTAT", ""), i.get("CLNDN", "")]
            w("\t".join(row) + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit("usage: vcf_to_table.py annotated.vcf.gz SAMPLE")
    main(sys.argv[1], sys.argv[2])
