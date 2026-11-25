#!/usr/bin/env bash
# Bioinformatic reproducibility: the same FASTQ processed twice from scratch must give
# identical variant records (headers are ignored because they carry dates and command lines).
set -euo pipefail
PROFILE=${1:-test,local}
for i in 1 2; do
  nextflow run main.nf -profile "$PROFILE" -ansi-log false --outdir "results_repro_$i" -w "work_repro_$i" > /dev/null
done
status=0
for f in variants/sampleA.g.vcf.gz variants/sampleA.filtered.vcf.gz annotated/sampleA.annotated.vcf.gz; do
  a=$(bcftools view -H "results_repro_1/sampleA/$f" | md5sum | cut -d' ' -f1)
  b=$(bcftools view -H "results_repro_2/sampleA/$f" | md5sum | cut -d' ' -f1)
  if [[ $a == "$b" ]]; then echo "identical  $f  $a"; else echo "DIFFERENT  $f"; status=1; fi
done
rm -rf results_repro_* work_repro_*
exit $status
