#!/usr/bin/env bash
# BQSR known-sites resources from the public GATK hg38 bundle, sliced to the panel.
set -euo pipefail
OUT=${1:-data}; PANEL=${PANEL:-assets/panel_regions.bed}
BUNDLE=https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0
mkdir -p "$OUT/known_sites"
bcftools view -R "$PANEL" -Oz -o "$OUT/known_sites/1000G_phase1.snps.panel.vcf.gz" \
  "$BUNDLE/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
bcftools view -R "$PANEL" -Oz -o "$OUT/known_sites/Mills_1000G.indels.panel.vcf.gz" \
  "$BUNDLE/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
for f in "$OUT"/known_sites/*.vcf.gz; do bcftools index -t -f "$f"; done
