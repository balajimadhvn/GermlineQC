# GermlineQC

[![CI](https://github.com/balajimadhvn/GermlineQC/actions/workflows/ci.yml/badge.svg)](https://github.com/balajimadhvn/GermlineQC/actions/workflows/ci.yml)
![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A524.04-23aa62)
![GATK](https://img.shields.io/badge/GATK-4.6.1.0-blue)
![License: MIT](https://img.shields.io/badge/license-MIT-lightgrey)

A germline short-variant pipeline written in Nextflow (FASTQ → BWA-MEM → GATK → annotated VCF) together
with the analytical validation a clinical lab needs before it can use one: accuracy against GIAB truth sets,
repeatability, between-run reproducibility and a depth-based limit of detection. Acceptance criteria were
set up front and the results are reported against them, including the one criterion the data did not meet.

Variants are then classified under the ACMG/AMP guidelines by a rule engine that uses the ClinGen points
system, and that engine is benchmarked against expert-reviewed ClinVar classifications. HLA typing of the
same GIAB sample lives in the companion repo [GIAB-HLA](https://github.com/balajimadhvn/GIAB-HLA).

```
FASTQ ─ fastp ─ BWA-MEM ─ MarkDuplicates ─ BQSR ─ HaplotypeCaller (GVCF) ─ GenotypeGVCFs
                                   │                                             │
                              mosdepth, samtools                       hard filters + bcftools norm
                                   │                                             │
                              QC gates (PASS/FAIL)                      SnpEff + ClinVar annotation
                                   └───────────────┬─────────────────────────────┘
                                         SQLite variant store, MultiQC
                                                   │
                         ACMG/AMP classification (gnomAD, REVEL, ClinGen, ClinVar) → benchmark vs ClinVar
                                                   │
                         validation: RTG vcfeval vs GIAB → Python metrics → R Markdown report
```

## Validation results

15 runs of GIAB reference material (HG001 and HG002, Illumina HiSeq 2500, 2 x 148 bp) over a 1.5 Mb chr20
panel covering the *JAG1* and *GNAS* loci, scored inside GIAB v4.2.1 high-confidence regions.


**[View the full validation report](https://balajimadhvn.github.io/GermlineQC/validation/validation_report.html)**
&nbsp;·&nbsp; [MultiQC report for all 15 runs](https://balajimadhvn.github.io/GermlineQC/validation/multiqc_report.html)
&nbsp;·&nbsp; [result tables](validation)

| Criterion | Observed | 95% CI lower bound | Target | Result |
|---|---|---|---|---|
| SNV sensitivity (8 runs at 30x, 14,493 truth SNVs) | 99.25% | 99.09% | ≥ 99% | PASS |
| SNV PPV | 99.89% | 99.82% | ≥ 99% | PASS |
| Indel sensitivity (2,605 truth indels) | 98.89% | 98.41% | ≥ 95% | PASS |
| Indel PPV | 99.63% | 99.32% | ≥ 95% | PASS |
| Specificity, per base (11.6 M reference bases) | 99.9998% | | ≥ 99.99% | PASS |
| Reproducibility, 4 runs on 2 instruments over 4 months (worst pair) | 98.01% | | ≥ 98% | PASS |
| Repeatability, 3 replicates from one batch (worst pair) | 97.92% | | ≥ 99% | **FAIL** |
| Bioinformatic reproducibility (same input, two clean runs) | identical VCF records | | identical | PASS |

<p align="center">
  <img src="validation/figures/titration-plot-1.png" width="49%" alt="Sensitivity by depth">
  <img src="validation/figures/accuracy-plot-1.png" width="49%" alt="Per-run accuracy">
</p>

**Why repeatability failed, and what changed.** Discordant calls between replicates are much more often
low-depth calls (22% have fewer than 15 reads, against 3% of concordant calls). The same 30x replicates also
miss the per-sample coverage gate: only about 90% of target bases reach 20x, against the required 95%. The
depth titration shows the gate is first met at 50x, where SNV and indel sensitivity are both 100% on this panel.
So the problem is read sampling at 30x, not the caller. The corrective action is to set the assay's
minimum depth to 50x (now enforced by the QC gate in `assets/qc_thresholds.yaml`) and repeat the precision
studies at that depth. The acceptance criteria were not
loosened after seeing the data. Details are in section 6 of the report.

## ACMG/AMP classification

`bin/acmg_classify.py` assigns each variant the ACMG/AMP criteria it can evaluate from public data, then
combines them with the ClinGen points system (Tavtigian et al. 2020: supporting 1, moderate 2, strong 4,
very strong 8; ≥ 10 pathogenic, 6 to 9 likely pathogenic, 0 to 5 VUS, −1 to −6 likely benign, ≤ −7 benign).

| Criterion | Evidence used | Rule |
|---|---|---|
| PVS1 | SnpEff consequence + ClinGen haploinsufficiency score | Null variant in a gene where loss of function causes disease; Strong in the last exon, Moderate for start-loss |
| PS1 / PM5 | ClinVar (≥ 1 star) | Same amino-acid change, or another missense change at the same residue, already classified P/LP |
| PM2 (Supporting) | gnomAD v4.1 exomes + genomes | Absent or below the inheritance-specific ceiling |
| PM4 | SnpEff | In-frame indel or stop-loss |
| PP3 / BP4 | REVEL | Pejaver et al. 2022 calibrated thresholds, Supporting to Strong |
| BA1 / BS1 | gnomAD filtering allele frequency (FAF95, ancestry max) | 5%, and a gene-inheritance-specific threshold |
| BP7 | SnpEff | Synonymous, outside splice regions |

All thresholds are in [`assets/acmg_rules.yaml`](assets/acmg_rules.yaml) and were fixed before the benchmark
was run. Criteria that need case-level information (segregation, de novo status, phenotype, functional
studies) are deliberately not scored, so the engine is a first-pass triage tool for a variant scientist and
not a sign-out.

### Benchmark against ClinVar

The truth set is every chr20 ClinVar germline SNV or small indel (≤ 50 bp) classified P/LP or B/LB with at
least 2 review stars: 9,471 variants (1,391 P/LP, 8,080 B/LB). Each variant's own ClinVar record is
excluded when PS1/PM5 are evaluated (`--exclude-self`). Without that, the truth label would leak into the
evidence.

| ClinVar | n | Same direction | VUS | Opposite direction |
|---|---|---|---|---|
| P/LP | 1,391 | 511 (36.7%, 34.2-39.3) | 878 (63.1%) | 2 (0.14%, 0.04-0.52) |
| B/LB | 8,080 | 6,098 (75.5%, 74.5-76.4) | 1,982 (24.5%) | 0 (0%, 0-0.05) |

95% Wilson intervals in brackets.

- **Safety:** 2 of 9,471 variants were called in the opposite direction to ClinVar, and no benign
  variant was called P/LP. Both errors are ClinVar P/LP missense variants (CFAP61, HNF4A) with low REVEL scores
  (0.12 and 0.16), so BP4 fired. These are the cases where a variant scientist's review catches what
  in-silico evidence misses.
- **Where it stops at VUS:** null variants reach LP in 45% of cases. Missense variants reach LP in 24%:
  without segregation, de novo or functional evidence, PM2 + PP3 reach at most 5 points. Intronic,
  splice-region and in-frame P/LP variants all stay VUS, because no splicing predictor is used.
- **Benign side:** synonymous (BP7) and missense (BP4, BS1) variants reach B/LB 81% and 85% of the time.
  The frequency criteria are not fully independent of ClinVar, whose submitters also read gnomAD, so read
  these as upper bounds.
- No variant reached Pathogenic. The Pathogenic tier needs very strong evidence plus case-level data that a
  rule engine cannot see.

**[View the ACMG benchmark report](https://balajimadhvn.github.io/GermlineQC/validation/acmg/acmg_report.html)**
&nbsp;·&nbsp; [result tables](validation/acmg)

## What is in the repo

| Path | What it does |
|---|---|
| [`main.nf`](main.nf), [`modules/`](modules) | Nextflow DSL2 pipeline (`GERMLINE` workflow plus a `VALIDATION` entry point) |
| [`conf/`](conf) | Resource limits, the bundled `test` profile and the GIAB `validation` profile |
| [`bin/qc_summary.py`](bin/qc_summary.py) | Collects fastp, samtools, Picard, mosdepth and bcftools metrics and applies the QC gates in [`assets/qc_thresholds.yaml`](assets/qc_thresholds.yaml) |
| [`bin/vcf_to_table.py`](bin/vcf_to_table.py) | Flattens the SnpEff + ClinVar VCF into a table (most severe consequence per variant) |
| [`bin/validation_metrics.py`](bin/validation_metrics.py) | Sensitivity, PPV, F1, per-base specificity, pairwise concordance, Wilson CIs, LoD, discordance analysis |
| [`bin/load_sqlite.py`](bin/load_sqlite.py), [`sql/`](sql) | SQLite schema, a `reportable` view and seven review queries ([results](validation/sql_query_results.txt)) |
| [`reports/validation_report.Rmd`](reports/validation_report.Rmd) | R Markdown validation report (ggplot2, dplyr, tidyr) |
| [`reports/cohort_summary.R`](reports/cohort_summary.R) | R script that queries the SQLite store through DBI/RSQLite |
| [`bin/acmg_classify.py`](bin/acmg_classify.py) | ACMG/AMP rule engine (ClinGen points), thresholds in [`assets/acmg_rules.yaml`](assets/acmg_rules.yaml) |
| [`bin/acmg_benchmark.py`](bin/acmg_benchmark.py), [`reports/acmg_report.Rmd`](reports/acmg_report.Rmd) | Concordance with ClinVar ≥ 2-star classifications, by consequence and by criterion |
| [`scripts/`](scripts) | GIAB data download, replicate design, determinism check, ACMG evidence download |
| [`docs/VALIDATION_PLAN.md`](docs/VALIDATION_PLAN.md) | Validation plan: scope, reference materials, metrics, acceptance criteria, sample size, revalidation triggers |
| [`docs/SOP.md`](docs/SOP.md) | Operating procedure for running the pipeline and reviewing QC before release |
| [`tests/`](tests) | pytest unit tests and an 80 kb test dataset of real HG001 reads |

## Quick start

```bash
git clone https://github.com/balajimadhvn/GermlineQC.git
cd GermlineQC
mamba env create -f environment.yml
conda activate germlineqc

# about 2 minutes: bundled 80 kb JAG1 slice, ~70x real HG001 reads
nextflow run main.nf -profile test,local
python -m pytest
```

Outputs land in `results_test/`: annotated VCF, BAM, per-sample QC table, `database/variants.db` and a
MultiQC report.

### Reproduce the validation

```bash
./scripts/fetch_giab.sh data            # panel slices of GIAB BAMs, truth sets, ClinVar, Ensembl GTF (~700 MB)
./scripts/fetch_known_sites.sh data     # BQSR resources from the GATK hg38 bundle
python scripts/make_validation_set.py --data data
nextflow run main.nf -entry VALIDATION -profile validation,local
```

The full run (15 samples) takes about an hour on a laptop with the default limit of 8 CPUs and 6.5 GB RAM
(`conf/base.config`). Raise the limits there on a bigger machine.

### Reproduce the ACMG benchmark

```bash
./scripts/fetch_acmg_resources.sh data  # ClinVar, gnomAD v4.1 and REVEL for chr20, ClinGen gene list
nextflow run main.nf -entry ACMG_BENCHMARK_RUN -profile acmg,local
```

Add `-profile acmg` to a normal run to classify every sample's variants; the classes are written to
`<sample>/acmg/<sample>.acmg.tsv` and loaded into the `acmg` table of the SQLite store.

### Your own data

```bash
# samplesheet.csv columns: sample,fastq_1,fastq_2
nextflow run main.nf -profile local \
  --input samplesheet.csv \
  --fasta genome.fa \
  --panel_bed targets.bed \
  --known_sites dbsnp.vcf.gz,mills.vcf.gz \
  --gtf Homo_sapiens.GRCh38.113.gtf \
  --clinvar clinvar.vcf.gz \
  --outdir results
```

## Design notes

* **Real replicates, not random splits.** GIAB's 300x HG001 data was sequenced on 12 flowcells across two
  HiSeq instruments between December 2013 and April 2014, and the flowcell ID survives in every read name.
  `make_validation_set.py` uses it to build reproducibility replicates from genuinely different sequencing runs.
  The script streams the FASTQs twice and never holds them in memory (peak RSS under 30 MB).
* **Haplotype-aware comparison.** Calls are scored with RTG vcfeval, the engine behind GA4GH benchmarking.
  An indel written two different ways still matches, and a call only counts as correct if the genotype is right.
* **Confidence intervals everywhere.** Proportions near 100% need Wilson score intervals. The validation plan
  includes the rule-of-three sample-size reasoning for how many variants are needed to claim 99%.
* **Hard filters instead of VQSR.** VQSR needs far more variants than a single-sample panel has, so SNVs and
  indels get GATK's recommended hard filters and multiallelic sites are split and left-aligned.
* **Annotation without a multi-GB database.** The SnpEff database is built inside the pipeline from the
  Ensembl GTF for the reference in use, and ClinVar significance, review status and disease are added with bcftools.
* **QC gates are data, not code.** Per-sample thresholds and validation criteria live in YAML files, so a
  threshold change is a reviewed diff with a history.

## Limitations

The panel is small (1.5 Mb), so indel estimates have wide intervals, and same-batch replicates top out at
about 34x on this data. The replicates are read subsets rather than separate library preparations, so wet-lab
variation is not measured. Copy-number, structural and mosaic variants are out of scope. These points are
discussed in the report and the validation plan.

## Data sources

GIAB HG001/HG002 reads and v4.2.1 benchmarks (NIST), GRCh38 chr20 (UCSC), Ensembl release 113 gene models,
ClinVar (NCBI), GATK hg38 resource bundle (Broad Institute), gnomAD v4.1 (Broad Institute), REVEL v1.3,
ClinGen gene curation list. All are public.

## License

MIT
