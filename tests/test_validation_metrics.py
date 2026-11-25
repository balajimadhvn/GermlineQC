import gzip
import math

import pytest

import validation_metrics as vm


def test_wilson_matches_reference_values():
    lo, hi = vm.wilson(95, 100)
    assert lo == pytest.approx(0.8882, abs=1e-4)
    assert hi == pytest.approx(0.9785, abs=1e-4)


def test_wilson_perfect_score_has_lower_bound_below_one():
    lo, hi = vm.wilson(300, 300)
    assert hi == pytest.approx(1.0) and 0.98 < lo < 1.0


def test_wilson_empty():
    assert all(math.isnan(x) for x in vm.wilson(0, 0))


def test_var_type():
    assert vm.var_type("A", "G") == "SNV"
    assert vm.var_type("A", "AT") == "INDEL"
    assert vm.var_type("A", "G,AT") == "INDEL"


def _vcf(path, records):
    with gzip.open(path, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n")
        for pos, ref, alt, gt in records:
            fh.write(f"chr20\t{pos}\t.\t{ref}\t{alt}\t50\tPASS\t.\tGT\t{gt}\n")


def test_accuracy_counts(tmp_path):
    d = tmp_path / "S1"
    d.mkdir()
    _vcf(d / "tp-baseline.vcf.gz", [(10, "A", "G", "0/1"), (20, "C", "T", "1/1"), (30, "AT", "A", "0/1")])
    _vcf(d / "tp.vcf.gz", [(10, "A", "G", "0/1"), (20, "C", "T", "1/1"), (30, "AT", "A", "0/1")])
    _vcf(d / "fn.vcf.gz", [(40, "G", "A", "0/1")])
    _vcf(d / "fp.vcf.gz", [(50, "T", "C", "0/1")])
    (d / "eval.bed").write_text("chr20\t0\t1000\n")
    rows = {r["type"]: r for r in vm.accuracy("S1", d)}
    assert rows["SNV"]["tp"] == 2 and rows["SNV"]["fn"] == 1 and rows["SNV"]["fp"] == 1
    assert rows["SNV"]["sensitivity"] == pytest.approx(2 / 3)
    assert rows["INDEL"]["sensitivity"] == 1.0
    # 1000 bases - 4 truth variants = 996 negatives, 1 FP
    assert rows["ALL"]["tn"] == 995
    assert rows["ALL"]["specificity"] == pytest.approx(995 / 996)


def test_concordance_requires_matching_genotype(tmp_path):
    a, b = tmp_path / "a.vcf.gz", tmp_path / "b.vcf.gz"
    _vcf(a, [(10, "A", "G", "0/1"), (20, "C", "T", "1/1"), (30, "G", "A", "0/1")])
    _vcf(b, [(10, "A", "G", "1/0"), (20, "C", "T", "0/1"), (40, "T", "C", "0/1")])
    regions = {"chr20": [(0, 1000)]}
    row = vm.concordance("repeatability", ["a", "b"], {"a": a, "b": b}, regions)[0]
    assert row["shared"] == 2 and row["genotype_match"] == 1
    assert row["concordance"] == pytest.approx(1 / 4)
