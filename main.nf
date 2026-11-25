#!/usr/bin/env nextflow
/*
 * GermlineQC
 * FASTQ -> fastp -> BWA-MEM -> MarkDuplicates -> BQSR -> HaplotypeCaller (GVCF)
 *       -> GenotypeGVCFs -> hard filters -> normalise -> SnpEff + ClinVar -> annotated VCF
 * plus coverage, QC acceptance checks, a SQLite variant store and MultiQC.
 *
 * The VALIDATION entry point benchmarks calls against GIAB truth sets.
 */
nextflow.enable.dsl = 2

include { PREPARE_REFERENCE   } from './modules/reference'
include { SNPEFF_BUILD        } from './modules/reference'
include { FASTP; FASTQC       } from './modules/qc'
include { BWA_MEM             } from './modules/align'
include { MARK_DUPLICATES     } from './modules/align'
include { BQSR                } from './modules/align'
include { ALIGN_STATS         } from './modules/align'
include { MOSDEPTH            } from './modules/align'
include { HAPLOTYPECALLER     } from './modules/call'
include { GENOTYPE_GVCF       } from './modules/call'
include { HARD_FILTER         } from './modules/call'
include { ANNOTATE            } from './modules/annotate'
include { QC_SUMMARY          } from './modules/report'
include { LOAD_SQLITE         } from './modules/report'
include { MULTIQC             } from './modules/report'
include { VCFEVAL             } from './modules/validate'
include { VALIDATION_METRICS  } from './modules/validate'
include { VALIDATION_REPORT   } from './modules/validate'

def helpMessage() {
    log.info """
    GermlineQC  v${workflow.manifest.version}

    Usage:
      nextflow run main.nf -profile test,local
      nextflow run main.nf --input samplesheet.csv --fasta ref.fa --panel_bed panel.bed -profile local
      nextflow run main.nf -entry VALIDATION -profile validation,local

    Key parameters:
      --input        CSV with columns sample,fastq_1,fastq_2
      --fasta        reference FASTA (indexes are built if missing)
      --panel_bed    target regions (BED)
      --known_sites  comma separated VCFs for BQSR (empty = skip BQSR)
      --gtf          Ensembl GTF gene models used to build the SnpEff database
      --clinvar      ClinVar VCF (bgzipped + tabix)
      --outdir       results directory [${params.outdir}]
    """.stripIndent()
}

def readSamplesheet(path) {
    Channel.fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            if (!row.sample || !row.fastq_1 || !row.fastq_2)
                error "Samplesheet rows need sample,fastq_1,fastq_2: ${row}"
            [ [id: row.sample], file(row.fastq_1, checkIfExists: true), file(row.fastq_2, checkIfExists: true) ]
        }
}

workflow GERMLINE {
    take:
    reads       // [meta, r1, r2]

    main:
    fasta     = file(params.fasta, checkIfExists: true)
    panel_bed = file(params.panel_bed, checkIfExists: true)

    PREPARE_REFERENCE(fasta)
    ref = PREPARE_REFERENCE.out.ref.collect()          // [fa, fai, dict, bwa index files...]

    FASTP(reads)
    FASTQC(reads)
    BWA_MEM(FASTP.out.reads, ref)
    MARK_DUPLICATES(BWA_MEM.out.bam)

    if (params.known_sites) {
        known = Channel.fromList(params.known_sites.tokenize(',').collect { [file(it, checkIfExists: true), file(it + '.tbi', checkIfExists: true)] })
                       .flatten().collect()
        BQSR(MARK_DUPLICATES.out.bam, ref, known, panel_bed)
        final_bam = BQSR.out.bam
    } else {
        final_bam = MARK_DUPLICATES.out.bam
    }

    ALIGN_STATS(final_bam, panel_bed)
    MOSDEPTH(final_bam, panel_bed)

    HAPLOTYPECALLER(final_bam, ref, panel_bed)
    GENOTYPE_GVCF(HAPLOTYPECALLER.out.gvcf, ref, panel_bed)
    HARD_FILTER(GENOTYPE_GVCF.out.vcf, ref)

    SNPEFF_BUILD(fasta, file(params.gtf, checkIfExists: true))
    ANNOTATE(HARD_FILTER.out.vcf, SNPEFF_BUILD.out.db, file(params.clinvar, checkIfExists: true),
             file(params.clinvar + '.tbi', checkIfExists: true))

    qc_inputs = FASTP.out.json
        .join(ALIGN_STATS.out.stats)
        .join(MARK_DUPLICATES.out.metrics)
        .join(MOSDEPTH.out.summary)
        .join(HARD_FILTER.out.stats)
    QC_SUMMARY(qc_inputs, file("${projectDir}/assets/qc_thresholds.yaml"))

    LOAD_SQLITE(ANNOTATE.out.tsv.map { it[1] }.collect(), QC_SUMMARY.out.json.map { it[1] }.collect())

    MULTIQC(
        FASTP.out.json.map { it[1] }
            .mix(FASTQC.out.zip.map { it[1] })
            .mix(ALIGN_STATS.out.all.map { it[1] })
            .mix(MARK_DUPLICATES.out.metrics.map { it[1] })
            .mix(MOSDEPTH.out.all.map { it[1] })
            .mix(HARD_FILTER.out.stats.map { it[1] })
            .mix(ANNOTATE.out.snpeff_csv.map { it[1] })
            .flatten().collect(),
        file("${projectDir}/assets/multiqc_config.yaml")
    )

    emit:
    vcf      = ANNOTATE.out.vcf
    filtered = HARD_FILTER.out.vcf
    qc       = QC_SUMMARY.out.json
}

workflow {
    if (params.help) { helpMessage(); exit 0 }
    if (!params.input) { helpMessage(); error "--input is required" }
    GERMLINE(readSamplesheet(params.input))
}

/*
 * Analytical validation: run the pipeline on GIAB reference material,
 * benchmark every sample against its truth set, then score accuracy,
 * repeatability, reproducibility and depth-based limit of detection.
 */
workflow VALIDATION {
    if (!params.validation_sheet) error "--validation_sheet is required"
    sheet = Channel.fromPath(params.validation_sheet, checkIfExists: true).splitCsv(header: true)

    reads = sheet.map { r -> [ [id: r.sample], file(r.fastq_1, checkIfExists: true), file(r.fastq_2, checkIfExists: true) ] }
    truth = sheet.map { r -> [ [id: r.sample], file(r.truth_vcf), file(r.truth_vcf + '.tbi'), file(r.truth_bed) ] }

    GERMLINE(reads)

    VCFEVAL(GERMLINE.out.filtered.join(truth), file(params.fasta), file(params.panel_bed))

    VALIDATION_METRICS(
        VCFEVAL.out.results.map { it[1] }.collect(),
        GERMLINE.out.qc.map { it[1] }.collect(),
        file(params.validation_sheet),
        file("${projectDir}/assets/acceptance_criteria.yaml")
    )
    VALIDATION_REPORT(VALIDATION_METRICS.out.tables, file("${projectDir}/reports/validation_report.Rmd"))
}
