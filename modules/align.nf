process BWA_MEM {
    tag "${meta.id}"
    label 'heavy'

    input:
    tuple val(meta), path(reads)
    path ref

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), emit: bam

    script:
    def rg = "@RG\\tID:${meta.id}\\tSM:${meta.id}\\tLB:${meta.id}_lib\\tPL:ILLUMINA\\tPU:${meta.id}.unit1"
    """
    bwa mem -t ${task.cpus} -K 100000000 -Y -R '${rg}' genome.fa ${reads[0]} ${reads[1]} \\
      | samtools sort -@ ${task.cpus} -m 512M -o ${meta.id}.sorted.bam -
    """
}

process MARK_DUPLICATES {
    tag "${meta.id}"
    label 'heavy'
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy', pattern: '*.metrics.txt'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.id}.md.bam"), path("${meta.id}.md.bai"), emit: bam
    tuple val(meta), path("${meta.id}.dup.metrics.txt"),                   emit: metrics

    script:
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.75) as int}m" MarkDuplicates \\
        -I ${bam} -O ${meta.id}.md.bam -M ${meta.id}.dup.metrics.txt \\
        --CREATE_INDEX true --VALIDATION_STRINGENCY SILENT --VERBOSITY ERROR
    """
}

process BQSR {
    tag "${meta.id}"
    label 'heavy'
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path known
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.recal.bam"), path("${meta.id}.recal.bai"), emit: bam
    path "${meta.id}.recal.table"

    script:
    def ks = known.findAll { it.name.endsWith('.vcf.gz') }.collect { "--known-sites ${it}" }.join(' ')
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.75) as int}m" BaseRecalibrator \\
        -R genome.fa -I ${bam} ${ks} -L ${panel_bed} -ip 100 -O ${meta.id}.recal.table
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.75) as int}m" ApplyBQSR \\
        -R genome.fa -I ${bam} --bqsr-recal-file ${meta.id}.recal.table \\
        -L ${panel_bed} -ip 100 -O ${meta.id}.recal.bam
    """
}

process ALIGN_STATS {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/alignment", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.flagstat.json"), path("${meta.id}.samtools_stats.txt"), emit: stats
    tuple val(meta), path("${meta.id}.{flagstat.json,samtools_stats.txt}"),                   emit: all

    script:
    """
    samtools flagstat -O json ${bam} > ${meta.id}.flagstat.json
    samtools stats -t ${panel_bed} ${bam} > ${meta.id}.samtools_stats.txt
    """
}

process MOSDEPTH {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/coverage", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.mosdepth.summary.txt"), path("${meta.id}.thresholds.bed.gz"), emit: summary
    tuple val(meta), path("${meta.id}.*"),                                                         emit: all

    script:
    """
    mosdepth --threads ${task.cpus} --no-per-base --by ${panel_bed} \\
        --thresholds 10,20,30,50 --mapq 20 ${meta.id} ${bam}
    """
}
