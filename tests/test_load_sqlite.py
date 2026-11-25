import json
import sqlite3
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_load_and_query(tmp_path):
    tsv = tmp_path / "S1.annotated.tsv"
    cols = ["sample", "chrom", "pos", "ref", "alt", "var_type", "filter", "qual", "genotype", "depth", "gq",
            "vaf", "gene", "effect", "impact", "transcript", "hgvs_c", "hgvs_p", "clinvar_id", "clinvar_sig",
            "clinvar_review", "clinvar_disease"]
    rows = [
        ["S1", "chr20", "100", "A", "G", "SNV", "PASS", "50", "het", "30", "99", "0.5", "JAG1",
         "missense_variant", "MODERATE", "ENST1", "c.1A>G", "p.Lys1Glu", "1", "Pathogenic", "", "Alagille"],
        ["S1", "chr20", "200", "C", "T", "SNV", "SNP_HARD", "10", "het", "5", "12", "0.2", "JAG1",
         "intron_variant", "MODIFIER", "ENST1", "", "", "", "", "", ""],
    ]
    tsv.write_text("\t".join(cols) + "\n" + "\n".join("\t".join(r) for r in rows) + "\n")
    qc = tmp_path / "S1.qc.json"
    qc.write_text(json.dumps({"sample": "S1", "overall": "PASS",
                              "metrics": [{"metric": "mean_target_depth", "value": 31.2, "min": 20,
                                           "max": None, "status": "PASS"}]}))
    db = tmp_path / "v.db"
    subprocess.run([sys.executable, ROOT / "bin" / "load_sqlite.py", "--schema", ROOT / "sql" / "schema.sql",
                    "--db", db, "--variants", tsv, "--qc", qc], check=True)

    con = sqlite3.connect(db)
    assert con.execute("SELECT COUNT(*) FROM variants").fetchone()[0] == 2
    reportable = con.execute("SELECT gene, hgvs_p FROM reportable").fetchall()
    assert reportable == [("JAG1", "p.Lys1Glu")]
    # every example query in sql/queries.sql must run against the schema
    con.executescript((ROOT / "sql" / "queries.sql").read_text())
