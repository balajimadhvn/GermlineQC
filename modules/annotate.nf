process ANNOTATE {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/annotated", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path snpeff_db
    path clinvar
    path clinvar_tbi

    output:
    tuple val(meta), path("${meta.id}.annotated.vcf.gz"), path("${meta.id}.annotated.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.annotated.tsv"),                                            emit: tsv
    tuple val(meta), path("${meta.id}.snpeff.csv"),                                               emit: snpeff_csv
    path "${meta.id}.snpeff.html"

    script:
    """
    snpEff -Xmx${(task.memory.toMega() * 0.75) as int}m ann -c ${snpeff_db}/snpEff.config -dataDir \$PWD/${snpeff_db}/data \\
        -csvStats ${meta.id}.snpeff.csv -stats ${meta.id}.snpeff.html \\
        -no-intergenic ${params.snpeff_genome} ${vcf} | bgzip > snpeff.vcf.gz
    tabix -p vcf snpeff.vcf.gz

    bcftools annotate -a ${clinvar} \\
        -c INFO/CLNSIG,INFO/CLNREVSTAT,INFO/CLNDN,INFO/CLNVC,INFO/CLNID:=ID \\
        -h <(printf '##INFO=<ID=CLNID,Number=1,Type=String,Description="ClinVar variation ID">\\n') \\
        snpeff.vcf.gz -Oz -o ${meta.id}.annotated.vcf.gz
    bcftools index -t ${meta.id}.annotated.vcf.gz

    vcf_to_table.py ${meta.id}.annotated.vcf.gz ${meta.id} > ${meta.id}.annotated.tsv
    """
}
