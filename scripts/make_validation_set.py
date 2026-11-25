#!/usr/bin/env python3
"""
Build the analytical-validation sample set from the GIAB 300x panel FASTQs.

The HG001 300x data set was sequenced on 12 HiSeq 2500 flowcells across two
instruments between Dec 2013 and Apr 2014. Read names keep the flowcell ID, so
replicates can be drawn from genuinely different sequencing runs:

  * repeatability   - 3 disjoint random subsets of one batch (same library, same instrument,
                      flowcells run within three weeks: Jan 15 - Feb 7 2014)
  * reproducibility - 4 replicates, each from a different run date / instrument
  * depth titration - pooled reads downsampled to 5x ... 100x (limit-of-detection by depth)
  * accuracy        - HG001 and HG002 at ~30x against their GIAB v4.2.1 truth sets

Every subset is deterministic (hash-based sampling with a fixed seed) and the
FASTQs are streamed twice, so memory use stays flat regardless of input size.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
from dataclasses import dataclass
from pathlib import Path

PANEL_BP = 1_500_000

# flowcells grouped into sequencing runs (date, instrument) -> flowcells on that run
RUNS = {
    "run_2013-12-19_D00360": ["H814YADXX", "H81VLADXX"],
    "run_2013-12-23_D00360": ["H88WKADXX", "H88U0ADXX"],
    "batch_2014-01_HISEQ1":  ["H8962ADXX", "H894YADXX", "H8GV6ADXX", "H8GVUADXX", "H8G92ADXX"],
    "run_2014-03-13_HISEQ1": ["H8GGVADXX", "H9258ADXX"],
    "run_2014-04-07_HISEQ1": ["H948VADXX", "H947YADXX"],
}
REPEATABILITY_RUN = "batch_2014-01_HISEQ1"
REPRODUCIBILITY_RUNS = ["run_2013-12-19_D00360", "run_2013-12-23_D00360",
                        "run_2014-03-13_HISEQ1", "run_2014-04-07_HISEQ1"]
TITRATION = [5, 10, 15, 20, 30, 50, 100]


def read_pairs(r1_path: Path, r2_path: Path):
    """Stream (name, r1_record, r2_record) without holding the file in memory."""
    with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
        while True:
            a = "".join(f1.readline() for _ in range(4))
            b = "".join(f2.readline() for _ in range(4))
            if not a:
                return
            name = a[1:a.index("\n")].split()[0].removesuffix("/1")
            yield name, a, b


def unit_hash(name: str, seed: str) -> float:
    """Map a read name to a reproducible value in [0, 1)."""
    h = hashlib.blake2b(f"{seed}:{name}".encode(), digest_size=8).digest()
    return int.from_bytes(h, "big") / 2**64


def pair_bases(r1: str, r2: str) -> int:
    return len(r1.split("\n")[1]) + len(r2.split("\n")[1])


@dataclass
class Subset:
    """One output sample: a membership rule (filled in after pass 1) and its writers."""
    sample: str
    truth: str
    study: str
    level: str
    source: str
    pool: str                      # which pool the reads come from
    target: float                  # target depth
    seed: str
    split: int | None = None       # for repeatability: which third (0-2) of the pool
    frac: float = 1.0
    pairs: int = 0
    bases: int = 0

    def wants(self, name: str) -> bool:
        if self.split is not None and int(unit_hash(name, "split") * 3) != self.split:
            return False
        return unit_hash(name, self.seed) < self.frac


def pools_for(name: str, sample_id: str) -> list[str]:
    pools = [f"{sample_id}:all"]
    if sample_id == "HG001":
        run = FC_TO_RUN.get(name.split(":")[2])
        if run:
            pools.append(f"HG001:{run}")
    return pools


FC_TO_RUN = {fc: r for r, fcs in RUNS.items() for fc in fcs}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--data", type=Path, default=Path("data"))
    ap.add_argument("--target-depth", type=float, default=30.0)
    ap.add_argument("--seed", default="germlineqc-v1")
    args = ap.parse_args()

    reads, truth = args.data / "reads", args.data / "truth"
    outdir = args.data / "validation_reads"
    outdir.mkdir(parents=True, exist_ok=True)
    inputs = {s: (reads / f"{s}_300x_R1.fastq.gz", reads / f"{s}_300x_R2.fastq.gz") for s in ("HG001", "HG002")}

    # ---- pass 1: bases per pool (and per third of the repeatability pool) ----
    pool_bases: dict[str, int] = {}
    for sid, (r1, r2) in inputs.items():
        for name, a, b in read_pairs(r1, r2):
            n = pair_bases(a, b)
            for p in pools_for(name, sid):
                pool_bases[p] = pool_bases.get(p, 0) + n
                if p == f"HG001:{REPEATABILITY_RUN}":
                    k = f"{p}:third{int(unit_hash(name, 'split') * 3)}"
                    pool_bases[k] = pool_bases.get(k, 0) + n
    for p, n in sorted(pool_bases.items()):
        print(f"  pool {p:<40} ~{n / PANEL_BP:6.1f}x")

    # ---- define subsets ----
    subsets: list[Subset] = []
    rep_pool = f"HG001:{REPEATABILITY_RUN}"
    for i in range(3):
        sub = Subset(f"HG001_repeat_{i + 1}", "HG001", "repeatability", f"rep{i + 1}", REPEATABILITY_RUN,
                     rep_pool, args.target_depth, args.seed, split=i)
        sub.frac = min(1.0, args.target_depth * PANEL_BP / pool_bases[f"{rep_pool}:third{i}"])
        subsets.append(sub)
    for i, run in enumerate(REPRODUCIBILITY_RUNS, 1):
        subsets.append(Subset(f"HG001_run_{i}", "HG001", "reproducibility", run, run,
                              f"HG001:{run}", args.target_depth, args.seed))
    for d in TITRATION:
        subsets.append(Subset(f"HG001_depth_{d:03d}x", "HG001", "titration", str(d), "all runs pooled",
                              "HG001:all", d, args.seed + f":d{d}"))
    subsets.append(Subset("HG002_acc", "HG002", "accuracy", "30", "NHGRI 300x",
                          "HG002:all", args.target_depth, args.seed))
    for sub in subsets:
        if sub.split is None:
            sub.frac = min(1.0, sub.target * PANEL_BP / pool_bases[sub.pool])

    # ---- pass 2: stream reads straight into every subset that wants them ----
    writers = {}
    for sub in subsets:
        writers[sub.sample] = (gzip.open(outdir / f"{sub.sample}_R1.fastq.gz", "wt", compresslevel=3),
                               gzip.open(outdir / f"{sub.sample}_R2.fastq.gz", "wt", compresslevel=3))
    try:
        for sid, (r1, r2) in inputs.items():
            mine = [s for s in subsets if s.truth == sid]
            for name, a, b in read_pairs(r1, r2):
                pools = pools_for(name, sid)
                for sub in mine:
                    if sub.pool in pools and sub.wants(name):
                        w1, w2 = writers[sub.sample]
                        w1.write(a)
                        w2.write(b)
                        sub.pairs += 1
                        sub.bases += pair_bases(a, b)
    finally:
        for w1, w2 in writers.values():
            w1.close()
            w2.close()

    rows = []
    for sub in subsets:
        depth = sub.bases / PANEL_BP
        print(f"  {sub.sample:<22} {sub.pairs:>9,} pairs  ~{depth:5.1f}x  {sub.source}")
        rows.append(dict(sample=sub.sample,
                         fastq_1=str((outdir / f"{sub.sample}_R1.fastq.gz").resolve()),
                         fastq_2=str((outdir / f"{sub.sample}_R2.fastq.gz").resolve()),
                         truth_vcf=str((truth / f"{sub.truth}.truth.vcf.gz").resolve()),
                         truth_bed=str((truth / f"{sub.truth}.confident.bed").resolve()),
                         reference_material=sub.truth, study=sub.study, level=sub.level,
                         nominal_depth=round(depth, 1), source=sub.source))

    sheet = args.data / "validation_samplesheet.csv"
    with sheet.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {sheet} ({len(rows)} samples)")


if __name__ == "__main__":
    main()
