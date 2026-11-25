#!/usr/bin/env bash
# Evidence sources for the ACMG/AMP classifier (chr20 only, to match the validation panel):
#   ClinVar          - lookups for PS1/PM5, and the 2-star+ benchmark set
#   gnomAD v4.1      - exomes + genomes population frequencies (BA1, BS1, PM2), sliced remotely
#   REVEL v1.3       - missense pathogenicity score (PP3 / BP4, ClinGen SVI calibrated thresholds)
#   ClinGen dosage   - haploinsufficiency scores (PVS1 applicability)
set -euo pipefail
OUT=${1:-data}/acmg
PANEL=${PANEL:-assets/panel_regions.bed}
mkdir -p "$OUT/tmp"

CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
GNOMAD=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf
REVEL=https://zenodo.org/records/7072866/files/revel-v1.3_all_chromosomes.zip
CLINGEN=https://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh38.tsv

echo "[clinvar] chr20"
bcftools view -r 20 "$CLINVAR" | bcftools annotate --rename-chrs <(echo -e "20\tchr20") -Oz -o "$OUT/clinvar.chr20.vcf.gz"
bcftools index -t -f "$OUT/clinvar.chr20.vcf.gz"

echo "[clinvar] benchmark set: >= 2 review stars, germline P/LP or B/LB, SNV or indel <= 50 bp"
# CLNREVSTAT holds commas, which bcftools expressions split into separate values, so filter the text with awk
bcftools view "$OUT/clinvar.chr20.vcf.gz" | awk -F'\t' '
  /^#/ {print; next}
  {
    rev = ""; sig = ""
    n = split($8, kv, ";")
    for (i = 1; i <= n; i++) {
      if (kv[i] ~ /^CLNREVSTAT=/) rev = substr(kv[i], 12)
      else if (kv[i] ~ /^CLNSIG=/) sig = substr(kv[i], 8)
    }
    stars2 = (rev == "criteria_provided,_multiple_submitters,_no_conflicts" || rev == "reviewed_by_expert_panel" || rev == "practice_guideline")
    known  = (sig ~ /^(Pathogenic|Likely_pathogenic|Pathogenic\/Likely_pathogenic|Benign|Likely_benign|Benign\/Likely_benign)$/)
    small  = (length($4) <= 51 && length($5) <= 51 && $5 != ".")
    if (stars2 && known && small) print
  }' | bgzip > "$OUT/clinvar_benchmark.vcf.gz"
bcftools index -t -f "$OUT/clinvar_benchmark.vcf.gz"
echo "  $(bcftools view -H "$OUT/clinvar_benchmark.vcf.gz" | wc -l) benchmark variants"

echo "[gnomAD] v4.1 exomes + genomes at panel regions and benchmark sites"
{ awk -v OFS='\t' '{print $1, $2+1, $3}' "$PANEL"
  bcftools query -f '%CHROM\t%POS\t%POS\n' "$OUT/clinvar_benchmark.vcf.gz"; } | sort -k1,1 -k2,2n -u > "$OUT/tmp/sites.tsv"
# one ranged request per site is slow over HTTP: read sites within 2 kb of each other as one window,
# then keep only the wanted positions with -T
awk -v OFS='\t' 'NR==1 {c=$1; s=$2; e=$3; next}
  $1!=c || $2-e > 2000 {print c, s, e; c=$1; s=$2}
  {if ($3 > e) e=$3}
  END {print c, s, e}' "$OUT/tmp/sites.tsv" > "$OUT/tmp/windows.tsv"
fetch_gnomad() {
  local src=$1
  [[ -s $OUT/gnomad.$src.vcf.gz.tbi ]] && return
  bcftools view -R "$OUT/tmp/windows.tsv" -T "$OUT/tmp/sites.tsv" "$GNOMAD/$src/gnomad.$src.v4.1.sites.chr20.vcf.bgz" \
    | bcftools annotate -x '^INFO/AC,^INFO/AN,^INFO/AF,^INFO/AF_grpmax,^INFO/fafmax_faf95_max,^INFO/nhomalt' \
    | bcftools view -i 'INFO/AC>0' -Oz -o "$OUT/gnomad.$src.vcf.gz"
  bcftools index -t -f "$OUT/gnomad.$src.vcf.gz"
  echo "  $src: $(bcftools view -H "$OUT/gnomad.$src.vcf.gz" | wc -l) records"
}
fetch_gnomad exomes & fetch_gnomad genomes & wait

echo "[REVEL] v1.3 chr20, GRCh38 coordinates"
if [[ ! -s $OUT/revel.chr20.tsv.gz ]]; then
  curl -fsSL -o "$OUT/tmp/revel.zip" "$REVEL"
  # one row per transcript in the source; keep the highest score per allele
  unzip -p "$OUT/tmp/revel.zip" \
    | awk -F, -v OFS='\t' 'NR==1 {next} $1=="20" && $3!="." {k = $3 OFS $4 OFS $5; if (!(k in s) || $8 > s[k]) s[k] = $8}
                           END {for (k in s) print "chr20", k, s[k]}' \
    | sort -k2,2n -k3,3 -k4,4 | bgzip > "$OUT/revel.chr20.tsv.gz"
  tabix -f -s1 -b2 -e2 "$OUT/revel.chr20.tsv.gz"
fi

echo "[ClinGen] gene curation list"
curl -fsSL "$CLINGEN" > "$OUT/clingen_genes.tsv"

rm -rf "$OUT/tmp"
echo "done -> $OUT"
