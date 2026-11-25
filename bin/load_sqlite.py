#!/usr/bin/env python3
"""Load annotated variant tables and QC results into a SQLite database."""
from __future__ import annotations

import argparse
import csv
import json
import sqlite3
from pathlib import Path

INT_COLS = {"pos", "depth", "gq"}
REAL_COLS = {"qual", "vaf"}


def cast(col: str, value: str):
    if value in ("", "."):
        return None
    if col in INT_COLS:
        return int(value)
    if col in REAL_COLS:
        return float(value)
    return value


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--schema", type=Path, required=True)
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--variants", type=Path, nargs="+", required=True)
    ap.add_argument("--qc", type=Path, nargs="+", required=True)
    a = ap.parse_args()

    con = sqlite3.connect(a.db)
    con.executescript(a.schema.read_text())

    for qc_file in a.qc:
        qc = json.loads(qc_file.read_text())
        con.execute("INSERT OR REPLACE INTO samples(sample_id, qc_status) VALUES (?, ?)",
                    (qc["sample"], qc["overall"]))
        con.executemany(
            "INSERT INTO qc_metrics(sample_id, metric, value, min_allowed, max_allowed, status) VALUES (?,?,?,?,?,?)",
            [(qc["sample"], m["metric"], m["value"], m["min"], m["max"], m["status"]) for m in qc["metrics"]])

    n = 0
    for tsv in a.variants:
        with tsv.open() as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            cols = reader.fieldnames
            sql = f"INSERT INTO variants({', '.join(c if c != 'sample' else 'sample_id' for c in cols)}) " \
                  f"VALUES ({', '.join('?' * len(cols))})"
            rows = [[cast(c, r[c]) for c in cols] for r in reader]
            con.executemany(sql, rows)
            n += len(rows)
    con.commit()
    print(f"loaded {n} variants from {len(a.variants)} samples into {a.db}")


if __name__ == "__main__":
    main()
