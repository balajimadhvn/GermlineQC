// SnpEff annotation for VCFs that do not come from the calling workflow (ClinVar lookups, benchmark set)
process SNPEFF_VCF {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(vcf)
    path snpeff_db

    output:
    tuple val(meta), path("${meta.id}.snpeff.vcf.gz"), emit: vcf

    script:
    """
    snpEff -Xmx${(task.memory.toMega() * 0.75) as int}m ann -c ${snpeff_db}/snpEff.config \\
        -dataDir \$PWD/${snpeff_db}/data -noStats -no-intergenic ${params.snpeff_genome} ${vcf} \\
        | bgzip > ${meta.id}.snpeff.vcf.gz
    """
}

process ACMG_CLASSIFY {
    tag "${meta.id}"
    label 'light'
    publishDir { meta.benchmark ? "${params.outdir}/acmg_benchmark" : "${params.outdir}/${meta.id}/acmg" }, mode: 'copy'

    input:
    tuple val(meta), path(vcf)
    path gnomad
    path revel
    path clingen
    path clinvar_annotated
    path rules

    output:
    tuple val(meta), path("${meta.id}.acmg.tsv"), emit: tsv

    script:
    def bench = meta.benchmark ? '--exclude-self --all-records' : ''
    def gn = gnomad.findAll { it.name.endsWith('.vcf.gz') }.join(' ')
    def rv = revel.find { it.name.endsWith('.tsv.gz') }
    """
    acmg_classify.py --vcf ${vcf} --sample ${meta.id} --gnomad ${gn} --revel ${rv} \\
        --clingen ${clingen} --clinvar-annotated ${clinvar_annotated} --rules ${rules} ${bench} \\
        --out ${meta.id}.acmg.tsv
    """
}

process ACMG_BENCHMARK {
    label 'light'
    publishDir "${params.outdir}/acmg_benchmark", mode: 'copy'

    input:
    path classified
    path truth_vcf

    output:
    path "tables", emit: tables

    script:
    """
    acmg_benchmark.py --classified ${classified} --truth ${truth_vcf} --outdir tables
    """
}

process ACMG_REPORT {
    label 'light'
    publishDir "${params.outdir}/acmg_benchmark", mode: 'copy'

    input:
    path tables
    path rmd

    output:
    path "acmg_report.html"
    path "figures"

    script:
    """
    cp -L ${rmd} report.Rmd
    Rscript -e "rmarkdown::render('report.Rmd', output_file = 'acmg_report.html', params = list(tables = '\$PWD/${tables}'), quiet = TRUE)"
    """
}
