#!/usr/bin/env bash
# Pull only the panel regions from public GIAB / NCBI / Ensembl / UCSC sources.
# Remote BAM/VCF slicing via htslib means we never download whole genomes.
# Total download is roughly 700 MB.
set -euo pipefail

OUT=${1:-data}
PANEL=${PANEL:-assets/panel_regions.bed}
mkdir -p "$OUT"/{ref,reads,truth,annotation,tmp}

GIAB=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab
HG001_BAM=$GIAB/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/NHGRI_Illumina300X_novoalign_bams/HG001.GRCh38_full_plus_hs38d1_analysis_set_minus_alts.300x.bam
HG002_BAM=$GIAB/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam
HG001_TRUTH=$GIAB/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark
HG002_TRUTH=$GIAB/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark
CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
GTF=https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
CHR20=https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz

REGIONS=$(awk '{printf "%s:%d-%d ", $1, $2+1, $3}' "$PANEL")

echo "[ref] chr20"
if [[ ! -s $OUT/ref/chr20.fa ]]; then
  curl -fsSL "$CHR20" | gunzip -c > "$OUT/ref/chr20.fa"
fi

echo "[truth] GIAB v4.2.1 benchmark calls + confident regions"
for s in HG001 HG002; do
  url_var=${s}_TRUTH; base=${!url_var}
  bed_url=$base.bed; [[ $s == HG002 ]] && bed_url=${base}_noinconsistent.bed
  bcftools view -R "$PANEL" -Oz -o "$OUT/truth/$s.truth.vcf.gz" "$base.vcf.gz"
  bcftools index -t "$OUT/truth/$s.truth.vcf.gz"
  # clip the confident-region BED to the panel
  curl -fsSL "$bed_url" | awk -v OFS='\t' 'NR==FNR{c[NR]=$1; s[NR]=$2; e[NR]=$3; n=NR; next}
      {for(i=1;i<=n;i++) if($1==c[i]){a=($2>s[i]?$2:s[i]); b=($3<e[i]?$3:e[i]); if(a<b) print $1,a,b}}' \
      "$PANEL" - | sort -k1,1 -k2,2n > "$OUT/truth/$s.confident.bed"
done

echo "[reads] slicing 300x GIAB BAMs to the panel and converting back to FASTQ"
for s in HG001 HG002; do
  url_var=${s}_BAM
  if [[ ! -s $OUT/reads/${s}_300x_R1.fastq.gz ]]; then
    samtools view -@4 -b -F 0x900 -o "$OUT/tmp/$s.slice.bam" "${!url_var}" $REGIONS
    samtools collate -@4 -u -O "$OUT/tmp/$s.slice.bam" "$OUT/tmp/$s.collate" \
      | samtools fastq -@4 -s /dev/null -0 /dev/null \
          -1 "$OUT/reads/${s}_300x_R1.fastq.gz" -2 "$OUT/reads/${s}_300x_R2.fastq.gz" -
    rm -f "$OUT/tmp/$s.slice.bam"
  fi
done

echo "[annotation] ClinVar (panel only) and Ensembl chr20 gene models"
bcftools view -r "$(awk '{printf "%s%s:%d-%d", (NR>1?",":""), substr($1,4), $2+1, $3}' "$PANEL")" "$CLINVAR" \
  | bcftools annotate --rename-chrs <(echo -e "20\tchr20") -Oz -o "$OUT/annotation/clinvar.panel.vcf.gz"
bcftools index -t "$OUT/annotation/clinvar.panel.vcf.gz"
curl -fsSL "$GTF" | gunzip -c | awk -F'\t' -v OFS='\t' '$1=="20" {$1="chr20"; print}' \
  > "$OUT/annotation/genes.chr20.gtf"

rm -rf "$OUT/tmp"
echo "done -> $OUT"
