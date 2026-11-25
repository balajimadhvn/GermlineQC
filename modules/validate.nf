// GA4GH-style small-variant benchmarking with RTG vcfeval (the engine behind hap.py --engine vcfeval).
process VCFEVAL {
    tag "${meta.id}"
    label 'heavy'
    publishDir "${params.outdir}/validation/vcfeval", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi), path(truth), path(truth_tbi), path(truth_bed)
    path fasta
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}"), emit: results

    script:
    """
    # evaluate only inside GIAB confident regions that fall within the panel
    sort -k1,1 -k2,2n ${panel_bed} > panel.sorted.bed
    awk -v OFS='\\t' 'NR==FNR{c[NR]=\$1;s[NR]=\$2;e[NR]=\$3;n=NR;next}
        {for(i=1;i<=n;i++) if(\$1==c[i]){a=(\$2>s[i]?\$2:s[i]);b=(\$3<e[i]?\$3:e[i]); if(a<b) print \$1,a,b}}' \\
        panel.sorted.bed ${truth_bed} > eval.bed

    rtg format -o ref.sdf ${fasta} > /dev/null
    rtg RTG_MEM=${(task.memory.toMega() * 0.75) as int}m vcfeval -b ${truth} -c ${vcf} -t ref.sdf -e eval.bed \\
        --ref-overlap --vcf-score-field QUAL -T ${task.cpus} -o ${meta.id} > /dev/null
    cp eval.bed ${meta.id}/eval.bed
    """
}

process VALIDATION_METRICS {
    label 'light'
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path vcfeval_dirs
    path qc_jsons
    path sheet
    path criteria

    output:
    path "tables", emit: tables

    script:
    """
    validation_metrics.py --vcfeval ${vcfeval_dirs} --sheet ${sheet} --criteria ${criteria} \\
        --qc ${qc_jsons} --outdir tables
    """
}

process VALIDATION_REPORT {
    label 'light'
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    path tables
    path rmd

    output:
    path "validation_report.html"
    path "figures"

    script:
    """
    cp -L ${rmd} report.Rmd
    Rscript -e "rmarkdown::render('report.Rmd', output_file = 'validation_report.html', params = list(tables = '\$PWD/${tables}'), quiet = TRUE)"
    """
}
