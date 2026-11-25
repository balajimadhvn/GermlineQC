# SOP: Running the germline variant pipeline

**Document ID:** GQC-SOP-001, version 1.0

## 1. Purpose

To produce an annotated VCF and a QC record for each sample in a germline panel sequencing run, in a way
that can be traced and repeated.

## 2. Before the run

1. Confirm the pipeline release is a tagged, validated version (`git describe --tags`).
2. Activate the locked environment: `conda activate germlineqc`.
3. Prepare `samplesheet.csv` with columns `sample,fastq_1,fastq_2`. Sample IDs must match the LIMS IDs.
4. Check the FASTQ checksums against the sequencer output (`md5sum -c`).

## 3. Running

```bash
nextflow run main.nf \
  --input samplesheet.csv \
  --fasta /ref/GRCh38/genome.fa \
  --panel_bed /ref/panel/panel_v1.bed \
  --known_sites /ref/bqsr/1000G.snps.vcf.gz,/ref/bqsr/Mills.indels.vcf.gz \
  --gtf /ref/annotation/Homo_sapiens.GRCh38.113.gtf \
  --clinvar /ref/annotation/clinvar.vcf.gz \
  --outdir runs/RUN_ID \
  -profile local
```

## 4. QC review (must be done before any result is released)

1. Open `runs/RUN_ID/multiqc/multiqc_report.html` and check for outlier samples.
2. For each sample, open `<sample>/qc/<sample>.qc.tsv`:
   * every `FAIL` means the sample is **not reportable**. Record the reason and request a re-run or re-sequencing;
   * every `WARN` needs a reviewer's comment in the run record.
3. Run query 2 in `sql/queries.sql` against `database/variants.db` to list all out-of-range metrics for the run.
4. Variants returned by query 6 (low GQ, low depth or skewed allele fraction) need IGV review and, if
   reportable, orthogonal confirmation (e.g. Sanger) before release.

## 5. Outputs to archive

| File | Retention |
|---|---|
| `<sample>/annotated/<sample>.annotated.vcf.gz` | permanent |
| `<sample>/alignment/<sample>.recal.bam` | as per lab policy |
| `<sample>/qc/` | permanent |
| `pipeline_info/` (trace, timeline, report) | permanent |
| `multiqc/` | permanent |

## 6. Deviations

Any change to parameters for a single run (for example a lower depth threshold) is a deviation. It must be
recorded with the reason and approved by the lab director before results are released.
