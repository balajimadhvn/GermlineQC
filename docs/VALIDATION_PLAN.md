# Analytical Validation Plan

**Assay:** Germline SNV and small indel detection (panel, Illumina short reads)
**Pipeline:** GermlineQC v1.0.0 (BWA-MEM 0.7.18, GATK 4.6.1.0, SnpEff 5.2, RTG Tools 3.12.1)
**Document ID:** GQC-VP-001, version 1.0

This plan follows the structure that CAP (Molecular Pathology checklist requirements for NGS test validation and
bioinformatics) and NABL (ISO 15189:2022, clause 7.3.3 on validation of examination methods) expect from a
laboratory developed NGS test. The pipeline-specific requirements follow the joint AMP/CAP recommendations for
validating NGS bioinformatics pipelines (Roy et al., *J Mol Diagn* 2018). The study design (accuracy,
precision, limit of detection, reference materials, sample size by the rule of three) borrows from the AMP/CAP
guidance on NGS panel validation (Jennings et al., *J Mol Diagn* 2017), applied here to germline variants.

## 1. Scope

| Item | In scope | Out of scope |
|---|---|---|
| Variant types | SNVs, insertions and deletions up to 50 bp | CNVs, SVs, repeat expansions, mosaic / somatic variants |
| Regions | Panel BED (`assets/panel_regions.bed`) inside GIAB high-confidence regions | Low-mappability and segmental duplications outside GIAB confident regions |
| Specimen | Germline DNA, Illumina paired-end reads | Long-read, FFPE-specific artefacts |
| Software | Pipeline steps from FASTQ to annotated VCF | Wet-lab steps, interpretation and reporting |

## 2. Reference materials

| Material | Source | Truth set |
|---|---|---|
| HG001 (NA12878) | GIAB, NIST RM 8398 | GIAB v4.2.1 GRCh38 benchmark VCF + BED |
| HG002 (NA24385) | GIAB, NIST RM 8391 | GIAB v4.2.1 GRCh38 benchmark VCF + BED (no-inconsistent) |

Using a second, unrelated genome checks that performance is not tuned to one sample.

## 3. Performance characteristics and acceptance criteria

All criteria are set **before** the data is analysed. The values live in `assets/acceptance_criteria.yaml`.

| Characteristic | Definition | Acceptance criterion |
|---|---|---|
| Analytical sensitivity (SNV) | TP / (TP + FN) | >= 99.0% |
| Analytical sensitivity (indel) | TP / (TP + FN) | >= 95.0% |
| Positive predictive value (SNV) | TP / (TP + FP) | >= 99.0% |
| Positive predictive value (indel) | TP / (TP + FP) | >= 95.0% |
| Specificity (per base) | TN / (TN + FP), over confident reference bases | >= 99.99% |
| Repeatability | Pairwise genotype concordance between replicates of the same library, instrument and batch | >= 99.0% |
| Reproducibility | Pairwise genotype concordance between independent runs (different instruments and dates) | >= 98.0% |
| Limit of detection | Lowest mean depth from which SNV sensitivity stays >= 99% | Reported; sets the minimum depth QC gate |
| Reportable range | Panel bases with depth >= 20x in every validation run | >= 95% of target bases |

A true positive needs the right allele **and** the right zygosity. Comparison is haplotype-aware (RTG vcfeval),
so the same variant written two different ways (e.g. a left- vs right-aligned indel) is not counted as an error.

Each estimate is reported with a 95% Wilson score confidence interval. The Wilson interval behaves well
near 100%, where the normal-approximation interval breaks down.

### Sample size

With *n* true variants and zero misses, the 95% lower confidence bound on sensitivity is about `1 - 3/n` (rule of
three). Showing sensitivity >= 99% at 95% confidence therefore needs at least about 300 true variants
per variant class. SNVs on this panel clear that comfortably. Indels do not, so the indel estimate is reported
with its interval and flagged as a limitation (see the report).

## 4. Study design

| Study | Design | Samples |
|---|---|---|
| Accuracy | Each 30x run vs GIAB truth | HG001 x 7 runs, HG002 x 1 |
| Repeatability (within-run) | Same library, same instrument, one 3-week batch. Three disjoint read subsets at 30x | HG001_repeat_1..3 |
| Reproducibility (between-run) | Four sequencing runs on two HiSeq instruments, Dec 2013 to Apr 2014, 30x each | HG001_run_1..4 |
| Depth titration / LoD | Pooled reads downsampled to 5, 10, 15, 20, 30, 50, 100x | HG001_depth_* |
| Bioinformatic reproducibility | Same FASTQ processed twice from scratch gives identical GVCF, filtered and annotated VCF records | `scripts/check_determinism.sh` |

Replicates are built by `scripts/make_validation_set.py` from flowcell IDs in the GIAB read names, so the
"different runs" really are different sequencing runs and not random splits of one run.

## 5. Quality metrics applied to every sample

See `assets/qc_thresholds.yaml`. A sample that fails any metric is not reported and the run is reviewed.

| Metric | Threshold | Why |
|---|---|---|
| Q30 bases | >= 80% | Sequencing quality |
| Reads kept after trimming | >= 90% | Adapter / quality problems |
| Mapped reads | >= 98% | Contamination, wrong species, wrong reference |
| Properly paired | >= 90% | Library insert problems |
| Duplication | <= 30% (warning) | Low library complexity |
| Mean target depth | >= 50x (raised from 20x after validation) | Sensitivity and precision (see LoD and report section 6) |
| Target bases >= 20x | >= 95% | Reportable range |
| Ts/Tv | 1.8 to 2.6 (warning) | Systematic false positives shift Ts/Tv |

## 6. Revalidation triggers

A change to any of the following needs at least a partial revalidation (accuracy on one reference
sample, plus a concordance check against the previous version):

* aligner, variant caller or annotation software version
* reference genome build or panel BED
* filtering thresholds
* sequencing chemistry or instrument

Each pipeline release is tagged in git. `results_validation/pipeline_info/` keeps the software versions,
parameters and execution trace for the audit trail.

## 7. Records

| Record | Location |
|---|---|
| Acceptance criteria | `assets/acceptance_criteria.yaml` |
| Sample design | `data/validation_samplesheet.csv` (generated) |
| Raw comparison output | `results_validation/validation/vcfeval/<sample>/` |
| Summary tables | `results_validation/validation/tables/` |
| Signed report | `results_validation/validation/validation_report.html` (copy in `validation/`) |
| ACMG benchmark | `results/acmg_benchmark/` (copy in `validation/acmg/`) |

## 8. Variant classification (ACMG/AMP)

`bin/acmg_classify.py` is validated separately from the calling pipeline, because its errors are different:
a correct genotype can still be given the wrong clinical class.

| Item | Choice |
|---|---|
| Truth | ClinVar records with at least two review stars (multiple submitters, no conflicts, or expert panel), germline, SNV or indel up to 50 bp, on chr20 |
| Collapsed classes | P/LP and B/LB on both sides; predicted VUS is counted separately |
| Independence | the variant's own ClinVar record is never used as evidence for itself (`--exclude-self`) |
| Metrics | same-direction rate, VUS rate and opposite-direction rate per truth class, with Wilson 95% intervals; the same broken down by consequence type and by criterion |
| Thresholds | `assets/acmg_rules.yaml`, taken from ClinGen SVI publications and not tuned to the benchmark |

The opposite-direction rate (a ClinVar P/LP variant called B/LB, or the reverse) is the safety metric and
every such call is listed for review in `discordant.tsv`. A VUS result is safe but unhelpful. The engine is a
triage step: every classification still goes to a variant scientist, who adds the case-level criteria the
engine does not score (PS2, PS3, PS4, PM3, PM6, PP1, PP4, BS2, BS3, BS4, BP2, BP5).

A change to any of the following needs the benchmark to be re-run and compared with the previous version:
a rule threshold, the gnomAD or REVEL release, the ClinGen dosage list, or the SnpEff transcript set.
