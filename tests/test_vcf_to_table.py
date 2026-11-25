import gzip

import vcf_to_table as v


def test_genotype_classes():
    assert v.genotype("0/1") == "het"
    assert v.genotype("1|0") == "het"
    assert v.genotype("1/1") == "hom_alt"
    assert v.genotype("0/0") == "hom_ref"
    assert v.genotype("./.") == "nocall"


def test_top_annotation_picks_most_severe():
    low = "A|synonymous_variant|LOW|GENE1|g1|transcript|T1|protein_coding|2/5|c.3A>G|p.Lys1Lys|||||"
    high = "A|stop_gained|HIGH|GENE1|g1|transcript|T2|protein_coding|2/5|c.4C>T|p.Gln2*|||||"
    assert v.top_annotation(f"{low},{high}")[1] == "stop_gained"


def test_main_writes_expected_columns(tmp_path, capsys):
    vcf = tmp_path / "x.vcf.gz"
    with gzip.open(vcf, "wt") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n")
        fh.write("chr20\t100\t.\tA\tG\t50\tPASS\t"
                 "ANN=G|missense_variant|MODERATE|JAG1|ENSG1|transcript|ENST1|protein_coding|3/26|c.10A>G|p.Lys4Glu|||||;"
                 "CLNSIG=Pathogenic;CLNID=12345\tGT:AD:DP:GQ\t0/1:10,10:20:99\n")
        fh.write("chr20\t200\t.\tAT\tA\t30\tINDEL_HARD\t.\tGT:AD:DP:GQ\t1/1:0,8:8:24\n")
    v.main(str(vcf), "S1")
    lines = capsys.readouterr().out.strip().split("\n")
    header = lines[0].split("\t")
    first = dict(zip(header, lines[1].split("\t")))
    second = dict(zip(header, lines[2].split("\t")))
    assert header == v.COLUMNS
    assert first["gene"] == "JAG1" and first["hgvs_p"] == "p.Lys4Glu"
    assert first["vaf"] == "0.500" and first["clinvar_sig"] == "Pathogenic"
    assert second["var_type"] == "INDEL" and second["genotype"] == "hom_alt"
