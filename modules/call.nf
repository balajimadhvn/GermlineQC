process HAPLOTYPECALLER {
    tag "${meta.id}"
    label 'heavy'
    publishDir "${params.outdir}/${meta.id}/variants", mode: 'copy'

    input:
    tuple val(meta), path(bam), path(bai)
    path ref
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.g.vcf.gz"), path("${meta.id}.g.vcf.gz.tbi"), emit: gvcf

    script:
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.75) as int}m" HaplotypeCaller \\
        -R genome.fa -I ${bam} -L ${panel_bed} -ip 50 \\
        -ERC GVCF --native-pair-hmm-threads ${task.cpus} \\
        -O ${meta.id}.g.vcf.gz
    """
}

process GENOTYPE_GVCF {
    tag "${meta.id}"
    label 'light'

    input:
    tuple val(meta), path(gvcf), path(tbi)
    path ref
    path panel_bed

    output:
    tuple val(meta), path("${meta.id}.raw.vcf.gz"), path("${meta.id}.raw.vcf.gz.tbi"), emit: vcf

    script:
    """
    gatk --java-options "-Xmx${(task.memory.toMega() * 0.75) as int}m" GenotypeGVCFs \\
        -R genome.fa -V ${gvcf} -L ${panel_bed} -O ${meta.id}.raw.vcf.gz
    """
}

// Single-sample panels are too small for VQSR, so apply GATK hard filters per variant type,
// then left-align / split multiallelics so downstream comparison is representation-independent.
process HARD_FILTER {
    tag "${meta.id}"
    label 'light'
    publishDir "${params.outdir}/${meta.id}/variants", mode: 'copy'

    input:
    tuple val(meta), path(vcf), path(tbi)
    path ref

    output:
    tuple val(meta), path("${meta.id}.filtered.vcf.gz"), path("${meta.id}.filtered.vcf.gz.tbi"), emit: vcf
    tuple val(meta), path("${meta.id}.bcftools_stats.txt"),                                     emit: stats

    script:
    """
    gatk SelectVariants -V ${vcf} -select-type SNP -O snps.vcf.gz
    gatk SelectVariants -V ${vcf} -select-type INDEL -select-type MIXED -O indels.vcf.gz
    gatk VariantFiltration -V snps.vcf.gz -filter "${params.snp_filter}" --filter-name SNP_HARD -O snps.f.vcf.gz
    gatk VariantFiltration -V indels.vcf.gz -filter "${params.indel_filter}" --filter-name INDEL_HARD -O indels.f.vcf.gz
    gatk MergeVcfs -I snps.f.vcf.gz -I indels.f.vcf.gz -O merged.vcf.gz
    bcftools norm -m -any -f genome.fa merged.vcf.gz -Oz -o ${meta.id}.filtered.vcf.gz
    bcftools index -t ${meta.id}.filtered.vcf.gz
    bcftools stats -f PASS -F genome.fa ${meta.id}.filtered.vcf.gz > ${meta.id}.bcftools_stats.txt
    """
}
