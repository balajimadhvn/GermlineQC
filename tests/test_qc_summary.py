import gzip

import qc_summary as q

CRITERIA = {
    "mean_target_depth": {"min": 20},
    "duplication_rate": {"max": 0.30, "level": "warn"},
    "ts_tv": {"min": 1.8, "max": 2.6},
}


def status(rows, name):
    return next(r["status"] for r in rows if r["metric"] == name)


def test_pass_fail_and_warn():
    rows = q.evaluate({"mean_target_depth": 35.0, "duplication_rate": 0.45, "ts_tv": 3.1, "snvs": 10}, CRITERIA)
    assert status(rows, "mean_target_depth") == "PASS"
    assert status(rows, "duplication_rate") == "WARN"
    assert status(rows, "ts_tv") == "FAIL"
    assert status(rows, "snvs") == "INFO"


def test_boundaries_are_inclusive():
    rows = q.evaluate({"mean_target_depth": 20.0, "ts_tv": 2.6}, CRITERIA)
    assert status(rows, "mean_target_depth") == "PASS"
    assert status(rows, "ts_tv") == "PASS"


def test_coverage_metrics(tmp_path):
    summary = tmp_path / "s.mosdepth.summary.txt"
    summary.write_text("chrom\tlength\tbases\tmean\tmin\tmax\nchr20\t100\t3000\t30\t0\t50\n"
                       "total_region\t200\t6000\t30.0\t0\t50\n")
    th = tmp_path / "s.thresholds.bed.gz"
    with gzip.open(th, "wt") as fh:
        fh.write("#chrom\tstart\tend\tregion\t10X\t20X\n")
        fh.write("chr20\t0\t100\tA\t100\t90\n")
        fh.write("chr20\t200\t300\tB\t100\t70\n")
    m = q.coverage_metrics(summary, th)
    assert m["mean_target_depth"] == 30.0
    assert m["pct_target_10x"] == 1.0
    assert m["pct_target_20x"] == 0.8


def test_variant_metrics(tmp_path):
    stats = tmp_path / "b.txt"
    stats.write_text("SN\t0\tnumber of SNPs:\t120\nSN\t0\tnumber of indels:\t15\n"
                     "TSTV\t0\t80\t40\t2.00\t80\t40\t2.00\n")
    assert q.variant_metrics(stats) == {"snvs": 120, "indels": 15, "ts_tv": 2.0}
