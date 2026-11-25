from pathlib import Path

import yaml

import acmg_classify as a

RULES = yaml.safe_load((Path(__file__).resolve().parents[1] / "assets" / "acmg_rules.yaml").read_text())
HI = {"JAG1": "3", "RECGENE": "30", "OTHER": "1"}


def var(effects, gene="JAG1", hgvs_p="", exon=(5, 26), pos=100, alt="T"):
    return a.Variant("chr20", pos, "C", alt, gene=gene, effects=set(effects), hgvs_p=hgvs_p, exon=exon)


def codes(met):
    return {c for c, _ in met}


def test_point_combination_thresholds():
    assert a.combine([("PVS1", "VeryStrong"), ("PM2", "Supporting")]) == (9, "Likely pathogenic")
    assert a.combine([("PVS1", "VeryStrong"), ("PS1", "Strong")]) == (12, "Pathogenic")
    assert a.combine([("PM2", "Supporting"), ("PP3", "Moderate")]) == (3, "Uncertain significance")
    assert a.combine([("BP4", "Supporting"), ("BP7", "Supporting")]) == (-2, "Likely benign")
    assert a.combine([("BS1", "Strong"), ("BP4", "Strong")]) == (-8, "Benign")
    assert a.combine([("BA1", "StandAlone")]) == (None, "Benign")


def test_novel_stop_gained_in_haploinsufficient_gene():
    met = a.evaluate(var({"stop_gained"}), RULES, {}, {}, HI, [], False)
    assert ("PVS1", "VeryStrong") in met and ("PM2", "Supporting") in met
    assert a.combine(met)[1] == "Likely pathogenic"


def test_pvs1_downgrades():
    last_exon = a.evaluate(var({"frameshift_variant"}, exon=(26, 26)), RULES, {}, {}, HI, [], False)
    assert ("PVS1", "Strong") in last_exon
    start = a.evaluate(var({"start_lost"}), RULES, {}, {}, HI, [], False)
    assert ("PVS1", "Moderate") in start
    # LoF is not an established mechanism for this gene
    assert "PVS1" not in codes(a.evaluate(var({"stop_gained"}, gene="OTHER"), RULES, {}, {}, HI, [], False))


def test_population_frequency():
    v = var({"missense_variant"})
    common = {v.key: dict(ac=5000, an=100000, af=0.06, af_grpmax=0.08, faf=0.07, nhomalt=100)}
    assert a.evaluate(v, RULES, common, {}, HI, [], False) == [("BA1", "StandAlone")]
    bs1 = {v.key: dict(ac=300, an=100000, af=0.003, af_grpmax=0.004, faf=0.0025, nhomalt=0)}
    met = a.evaluate(v, RULES, bs1, {}, HI, [], False)
    assert "BS1" in codes(met) and "PM2" not in codes(met)
    # same frequency in a recessive gene is below the carrier-frequency threshold
    vr = var({"missense_variant"}, gene="RECGENE")
    assert "BS1" not in codes(a.evaluate(vr, RULES, {vr.key: bs1[v.key]}, {}, HI, [], False))


def test_revel_calibrated_strengths():
    v = var({"missense_variant"})
    assert ("PP3", "Strong") in a.evaluate(v, RULES, {}, {v.key: 0.95}, HI, [], False)
    assert ("PP3", "Supporting") in a.evaluate(v, RULES, {}, {v.key: 0.70}, HI, [], False)
    assert ("BP4", "Moderate") in a.evaluate(v, RULES, {}, {v.key: 0.10}, HI, [], False)
    assert not {"PP3", "BP4"} & codes(a.evaluate(v, RULES, {}, {v.key: 0.5}, HI, [], False))


def test_ps1_and_pm5_use_other_clinvar_variants_only():
    v = var({"missense_variant"}, hgvs_p="p.Arg184His")
    same_aa_other_nt = [(("chr20", 101, "G", "A"), ("JAG1", "Arg184"), "p.Arg184His")]
    other_aa = [(("chr20", 100, "C", "G"), ("JAG1", "Arg184"), "p.Arg184Pro")]
    itself = [(v.key, ("JAG1", "Arg184"), "p.Arg184His")]
    assert "PS1" in codes(a.evaluate(v, RULES, {}, {}, HI, same_aa_other_nt, False))
    assert "PM5" in codes(a.evaluate(v, RULES, {}, {}, HI, other_aa, False))
    assert not {"PS1", "PM5"} & codes(a.evaluate(v, RULES, {}, {}, HI, itself, True))


def test_inframe_and_synonymous():
    assert "PM4" in codes(a.evaluate(var({"inframe_deletion"}), RULES, {}, {}, HI, [], False))
    assert "BP7" in codes(a.evaluate(var({"synonymous_variant"}), RULES, {}, {}, HI, [], False))
    assert "BP7" not in codes(a.evaluate(var({"synonymous_variant", "splice_region_variant"}), RULES, {}, {}, HI, [], False))


def test_best_annotation_prefers_protein_coding_and_severity():
    ann = ("T|intron_variant|MODIFIER|JAG1|g|transcript|ENST1|retained_intron|||||||,"
           "T|missense_variant|MODERATE|JAG1|g|transcript|ENST2|protein_coding|3/26|c.10C>T|p.Arg4Cys|||||,"
           "T|stop_gained|HIGH|JAG1|g|transcript|ENST3|protein_coding|4/26|c.20C>T|p.Arg7*|||||")
    best = a.best_annotation("T", ann)
    assert best["transcript"] == "ENST3" and "stop_gained" in best["effects"]
