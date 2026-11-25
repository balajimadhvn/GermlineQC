process FASTP {
    tag "${meta.id}"
    label 'heavy'
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy', pattern: '*.{json,html}'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("${meta.id}_R{1,2}.trim.fastq.gz"), emit: reads
    tuple val(meta), path("${meta.id}.fastp.json"),           emit: json
    path "${meta.id}.fastp.html"

    script:
    """
    fastp -i ${r1} -I ${r2} \\
        -o ${meta.id}_R1.trim.fastq.gz -O ${meta.id}_R2.trim.fastq.gz \\
        --detect_adapter_for_pe --cut_tail --length_required 36 \\
        --thread ${task.cpus} \\
        --json ${meta.id}.fastp.json --html ${meta.id}.fastp.html \\
        --report_title "${meta.id}" 2> ${meta.id}.fastp.log
    """
}

process FASTQC {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/qc", mode: 'copy'

    input:
    tuple val(meta), path(r1), path(r2)

    output:
    tuple val(meta), path("*_fastqc.zip"), emit: zip
    path "*_fastqc.html"

    script:
    """
    fastqc --quiet --threads ${task.cpus} ${r1} ${r2}
    """
}
