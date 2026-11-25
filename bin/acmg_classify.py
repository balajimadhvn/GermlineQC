#!/usr/bin/env python3
"""
Rule-based ACMG/AMP germline variant classification (Richards et al. 2015) for the criteria that can be
evaluated from population data, in-silico predictors, variant consequence and ClinVar. Evidence is combined
with the ClinGen SVI point system (Tavtigian et al. 2020):

    supporting = 1, moderate = 2, strong = 4, very strong = 8 (benign criteria score negative)
    >= 10 Pathogenic | 6..9 Likely pathogenic | 0..5 VUS | -1..-6 Likely benign | <= -7 Benign
    BA1 is stand-alone Benign.

Criteria implemented (thresholds in assets/acmg_rules.yaml):

    PVS1   null variant in a gene where loss of function is a known mechanism (ClinGen HI score 3 or 30);
           strength reduced for last-exon truncations (NMD escape) and start-loss (Abou Tayoun et al. 2018, simplified)
    PS1    same amino-acid change as a ClinVar pathogenic missense variant
    PM5    different missense change at a residue with a ClinVar pathogenic missense variant
    PM2    absent or extremely rare in gnomAD v4.1 (applied at supporting strength, ClinGen SVI 2020)
    PM4    in-frame indel or stop-loss changing protein length
    PP3    REVEL above the ClinGen-calibrated thresholds (Pejaver et al. 2022)
    BA1    gnomAD filtering allele frequency >= 5%
    BS1    gnomAD filtering allele frequency above the expected disease frequency (inheritance-dependent)
    BP4    REVEL below the calibrated benign thresholds
    BP7    synonymous variant outside the splice region

Not implemented (need curated or case-level data): PS2, PS3, PS4, PM1, PM3, PM6, PP1, PP2, PP4, BS2, BS3,
BS4, BP1, BP2, BP3, BP5. Results are a research aid, not a clinical report.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

import yaml

STRENGTH_POINTS = {"Supporting": 1, "Moderate": 2, "Strong": 4, "VeryStrong": 8}
NULL_EFFECTS = {"stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant", "start_lost"}
INFRAME_EFFECTS = {"inframe_insertion", "inframe_deletion", "disruptive_inframe_insertion",
                   "disruptive_inframe_deletion", "stop_lost"}
SEVERITY = ["stop_gained", "frameshift_variant", "splice_acceptor_variant", "splice_donor_variant",
            "start_lost", "stop_lost", "disruptive_inframe_deletion", "disruptive_inframe_insertion",
            "inframe_deletion", "inframe_insertion", "missense_variant", "splice_region_variant",
            "synonymous_variant", "5_prime_UTR_variant", "3_prime_UTR_variant", "intron_variant"]
AA3 = r"[A-Z][a-z]{2}"


@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str = ""
    effects: set[str] = field(default_factory=set)
    transcript: str = ""
    hgvs_c: str = ""
    hgvs_p: str = ""
    exon: tuple[int, int] | None = None
    sample_gt: str = ""

    @property
    def key(self):
        return (self.chrom, self.pos, self.ref, self.alt)

    @property
    def residue(self) -> tuple[str, str] | None:
        """(gene, 'Arg123') for a missense change, used for PS1/PM5 lookups."""
        m = re.match(rf"p\.({AA3}\d+)({AA3})$", self.hgvs_p)
        return (self.gene, m.group(1)) if m else None


# ---------------------------------------------------------------- input parsing

def severity_rank(effects: set[str]) -> int:
    ranks = [SEVERITY.index(e) for e in effects if e in SEVERITY]
    return min(ranks) if ranks else len(SEVERITY)


def best_annotation(alt: str, ann: str) -> dict | None:
    """Most severe SnpEff annotation on a protein-coding transcript for this ALT allele."""
    best = None
    for entry in ann.split(","):
        f = entry.split("|")
        if len(f) < 11 or f[0] != alt:
            continue
        effects = set(f[1].split("&"))
        cand = dict(effects=effects, gene=f[3], biotype=f[7], transcript=f[6], rank=f[8],
                    hgvs_c=f[9], hgvs_p=f[10], sev=severity_rank(effects))
        if best is None or (cand["biotype"] == "protein_coding", -cand["sev"]) > \
                (best["biotype"] == "protein_coding", -best["sev"]):
            best = cand
    return best


def read_variants(vcf: Path, pass_only: bool) -> list[Variant]:
    out = []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if pass_only and f[6] not in ("PASS", "."):
                continue
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            for alt in f[4].split(","):
                if alt in ("*", "."):
                    continue
                v = Variant(f[0], int(f[1]), f[3], alt)
                if len(f) > 9:
                    v.sample_gt = f[9].split(":")[0]
                    if v.sample_gt.replace("|", "/") in ("0/0", "./.", "."):
                        continue
                a = best_annotation(alt, info.get("ANN", ""))
                if a:
                    v.gene, v.effects, v.transcript = a["gene"], a["effects"], a["transcript"]
                    v.hgvs_c, v.hgvs_p = a["hgvs_c"], a["hgvs_p"]
                    if "/" in a["rank"]:
                        x, y = a["rank"].split("/")
                        v.exon = (int(x), int(y))
                out.append(v)
    return out


def read_gnomad(paths: list[Path]) -> dict[tuple, dict]:
    """Highest-evidence frequency record per allele across gnomAD exomes and genomes."""
    freq: dict[tuple, dict] = {}
    for p in paths:
        with gzip.open(p, "rt") as fh:
            for line in fh:
                if line[0] == "#":
                    continue
                f = line.rstrip("\n").split("\t")
                info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
                for i, alt in enumerate(f[4].split(",")):
                    def val(k, cast=float):
                        x = info.get(k, ".").split(",")
                        return cast(x[i]) if i < len(x) and x[i] not in (".", "") else None
                    rec = dict(ac=val("AC", int) or 0, an=int(info.get("AN", 0) or 0), af=val("AF") or 0.0,
                               af_grpmax=val("AF_grpmax"), faf=val("fafmax_faf95_max"), nhomalt=val("nhomalt", int))
                    k = (f[0], int(f[1]), f[3], alt)
                    old = freq.get(k)
                    # keep the dataset with more called alleles (usually exomes for coding sites)
                    if old is None or rec["an"] > old["an"]:
                        freq[k] = rec
    return freq


def read_revel(path: Path, variants: list[Variant]) -> dict[tuple, float]:
    sites = sorted({(v.chrom, v.pos) for v in variants if "missense_variant" in v.effects})
    if not sites:
        return {}
    regions = "".join(f"{c}\t{p}\t{p}\n" for c, p in sites)
    res = subprocess.run(["tabix", "-R", "/dev/stdin", str(path)], input=regions, capture_output=True,
                         text=True, check=True)
    scores: dict[tuple, float] = {}
    for line in res.stdout.splitlines():
        c, p, ref, alt, score = line.split("\t")[:5]
        k = (c, int(p), ref, alt)
        scores[k] = max(scores.get(k, 0.0), float(score))
    return scores


def read_clingen(path: Path) -> dict[str, str]:
    """Gene symbol -> haploinsufficiency score ('3', '30', '2', ...)."""
    hi = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) > 4:
                hi[f[0]] = f[4]
    return hi


PATHOGENIC = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}


def read_clinvar_pathogenic_missense(annotated_clinvar: Path) -> list[tuple[tuple, str, str]]:
    """[(variant key, gene/residue, hgvs_p)] for ClinVar P/LP missense variants with >= 1 review star."""
    out = []
    for v, info in _iter_clinvar(annotated_clinvar):
        if info.get("CLNSIG") in PATHOGENIC and "no_assertion" not in info.get("CLNREVSTAT", "no_assertion") \
                and "missense_variant" in v.effects and v.residue:
            out.append((v.key, v.residue, v.hgvs_p))
    return out


def _iter_clinvar(path: Path):
    with gzip.open(path, "rt") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            info = dict(x.split("=", 1) for x in f[7].split(";") if "=" in x)
            a = best_annotation(f[4], info.get("ANN", ""))
            v = Variant(f[0], int(f[1]), f[3], f[4])
            if a:
                v.gene, v.effects, v.hgvs_p, v.hgvs_c = a["gene"], a["effects"], a["hgvs_p"], a["hgvs_c"]
            yield v, info


# ---------------------------------------------------------------- criteria

def inheritance(hi_score: str | None) -> str:
    if hi_score == "3":
        return "dominant"
    if hi_score == "30":
        return "recessive"
    return "unknown"


def evaluate(v: Variant, rules: dict, gnomad: dict, revel: dict, hi: dict,
             clinvar_missense: list, exclude_self: bool) -> list[tuple[str, str]]:
    """Return triggered criteria as (code, strength)."""
    met: list[tuple[str, str]] = []
    g = gnomad.get(v.key)
    faf = None
    if g:
        faf = g["faf"] if g["faf"] is not None else g["af_grpmax"] if g["af_grpmax"] is not None else g["af"]
    mode = inheritance(hi.get(v.gene))

    # ---- population data
    if faf is not None and faf >= rules["BA1_faf"]:
        return [("BA1", "StandAlone")]
    if faf is not None and faf >= rules["BS1_faf"][mode]:
        met.append(("BS1", "Strong"))
    if g is None or g["ac"] == 0 or (g["af_grpmax"] or g["af"]) <= rules["PM2_max_af"][mode]:
        met.append(("PM2", "Supporting"))

    # ---- null variants
    null = v.effects & NULL_EFFECTS
    if null and hi.get(v.gene) in ("3", "30"):
        if null == {"start_lost"}:
            met.append(("PVS1", "Moderate"))
        elif v.exon and v.exon[0] == v.exon[1] and null & {"stop_gained", "frameshift_variant"}:
            met.append(("PVS1", "Strong"))     # last exon: escapes NMD
        else:
            met.append(("PVS1", "VeryStrong"))

    # ---- protein length change
    if v.effects & INFRAME_EFFECTS and not null:
        met.append(("PM4", "Moderate"))

    # ---- missense: in-silico score and known pathogenic residues
    if "missense_variant" in v.effects and not null:
        score = revel.get(v.key)
        if score is not None:
            for strength, cut in rules["PP3_revel"].items():
                if score >= cut:
                    met.append(("PP3", strength))
                    break
            else:
                for strength, cut in rules["BP4_revel"].items():
                    if score <= cut:
                        met.append(("BP4", strength))
                        break
        if v.residue:
            others = [c for c in clinvar_missense if not (exclude_self and c[0] == v.key)]
            if any(r == v.residue and p == v.hgvs_p and k != v.key for k, r, p in others):
                met.append(("PS1", "Strong"))
            elif any(r == v.residue and p != v.hgvs_p for k, r, p in others):
                met.append(("PM5", "Moderate"))

    # ---- synonymous
    if v.effects and v.effects <= {"synonymous_variant"}:
        met.append(("BP7", "Supporting"))
    return met


def combine(met: list[tuple[str, str]]) -> tuple[int | None, str]:
    if any(c == "BA1" for c, _ in met):
        return None, "Benign"
    pts = 0
    for code, strength in met:
        p = STRENGTH_POINTS[strength]
        pts += -p if code.startswith("B") else p
    if pts >= 10:
        cls = "Pathogenic"
    elif pts >= 6:
        cls = "Likely pathogenic"
    elif pts >= 0:
        cls = "Uncertain significance"
    elif pts >= -6:
        cls = "Likely benign"
    else:
        cls = "Benign"
    return pts, cls


def fmt_criteria(met: list[tuple[str, str]]) -> str:
    default = {"PVS1": "VeryStrong", "PS1": "Strong", "BS1": "Strong", "PM4": "Moderate", "PM5": "Moderate",
               "BP7": "Supporting", "BA1": "StandAlone"}
    return ",".join(c if default.get(c) == s else f"{c}_{s}" for c, s in met)


# ---------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--vcf", type=Path, required=True, help="SnpEff-annotated VCF")
    ap.add_argument("--sample", required=True)
    ap.add_argument("--gnomad", type=Path, nargs="+", required=True)
    ap.add_argument("--revel", type=Path, required=True, help="bgzipped + tabix indexed chrom,pos,ref,alt,score")
    ap.add_argument("--clingen", type=Path, required=True)
    ap.add_argument("--clinvar-annotated", type=Path, required=True, help="SnpEff-annotated ClinVar VCF")
    ap.add_argument("--rules", type=Path, required=True)
    ap.add_argument("--exclude-self", action="store_true",
                    help="benchmark mode: never use a variant's own ClinVar record as evidence")
    ap.add_argument("--all-records", action="store_true", help="also classify non-PASS records")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    rules = yaml.safe_load(a.rules.read_text())
    variants = read_variants(a.vcf, pass_only=not a.all_records)
    gnomad = read_gnomad(a.gnomad)
    revel = read_revel(a.revel, variants)
    hi = read_clingen(a.clingen)
    clinvar_missense = read_clinvar_pathogenic_missense(a.clinvar_annotated)

    cols = ["sample", "chrom", "pos", "ref", "alt", "genotype", "gene", "transcript", "consequence", "hgvs_c",
            "hgvs_p", "gnomad_faf", "gnomad_af", "revel", "clingen_hi", "criteria", "points", "acmg_class"]
    counts: dict[str, int] = {}
    with a.out.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(cols)
        for v in variants:
            met = evaluate(v, rules, gnomad, revel, hi, clinvar_missense, a.exclude_self)
            pts, cls = combine(met)
            counts[cls] = counts.get(cls, 0) + 1
            g = gnomad.get(v.key) or {}
            w.writerow([a.sample, v.chrom, v.pos, v.ref, v.alt, v.sample_gt, v.gene, v.transcript,
                        "&".join(sorted(v.effects, key=lambda e: SEVERITY.index(e) if e in SEVERITY else 99)),
                        v.hgvs_c, v.hgvs_p,
                        "" if g.get("faf") is None else f"{g['faf']:.6g}",
                        "" if not g else f"{g['af']:.6g}",
                        "" if v.key not in revel else f"{revel[v.key]:.3f}",
                        hi.get(v.gene, ""), fmt_criteria(met), "" if pts is None else pts, cls])
    summary = ", ".join(f"{k}: {n}" for k, n in sorted(counts.items()))
    print(f"{a.sample}: {len(variants)} variants classified ({summary})", file=sys.stderr)


if __name__ == "__main__":
    main()
