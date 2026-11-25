process QC_SUMMARY {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy'

    input:
    tuple val(meta), path(fastp_json), path(flagstat), path(samstats), path(dup_metrics),
          path(mosdepth_summary), path(thresholds), path(bcfstats)
    path qc_thresholds

    output:
    tuple val(meta), path("${meta.id}.qc.json"), emit: json
    path "${meta.id}.qc.tsv"

    script:
    """
    qc_summary.py --sample ${meta.id} --fastp ${fastp_json} --flagstat ${flagstat} \\
        --dup ${dup_metrics} --mosdepth ${mosdepth_summary} --thresholds ${thresholds} \\
        --bcfstats ${bcfstats} --criteria ${qc_thresholds} \\
        --json ${meta.id}.qc.json --tsv ${meta.id}.qc.tsv
    """
}

process LOAD_SQLITE {
    label 'light'
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    path tsvs
    path qc_jsons
    path acmg

    output:
    path "variants.db",          emit: db
    path "example_queries.txt"

    script:
    """
    acmg_opt=\$(ls *.acmg.tsv 2>/dev/null | xargs -r echo --acmg)
    load_sqlite.py --schema ${projectDir}/sql/schema.sql --db variants.db --variants ${tsvs} --qc ${qc_jsons} \$acmg_opt
    sqlite3 -header -column variants.db < ${projectDir}/sql/queries.sql > example_queries.txt
    """
}

process MULTIQC {
    label 'light'
    publishDir "${params.outdir}/multiqc", mode: 'copy'

    input:
    path files
    path config

    output:
    path "multiqc_report.html"
    path "multiqc_report_data"

    script:
    """
    multiqc --force --config ${config} --filename multiqc_report.html . > multiqc.log 2>&1
    """
}
