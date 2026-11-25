process PREPARE_REFERENCE {
    tag "${fasta.simpleName}"
    label 'light'
    storeDir "${params.outdir}/reference"

    input:
    path fasta

    output:
    path "ref/*", emit: ref

    script:
    """
    mkdir -p ref
    cp -L ${fasta} ref/genome.fa
    samtools faidx ref/genome.fa
    gatk CreateSequenceDictionary -R ref/genome.fa -O ref/genome.dict --VERBOSITY ERROR
    bwa index ref/genome.fa 2> ref/bwa_index.log
    """
}

// Build a SnpEff database from Ensembl GTF gene models so annotation needs no multi-GB download.
process SNPEFF_BUILD {
    tag "${params.snpeff_genome}"
    label 'light'

    input:
    path fasta
    path gtf

    output:
    path "snpeff_db", emit: db

    script:
    def g = params.snpeff_genome
    """
    mkdir -p snpeff_db/data/${g}
    cp -L ${fasta} snpeff_db/data/${g}/sequences.fa
    cp -L ${gtf} snpeff_db/data/${g}/genes.gtf
    printf "data.dir = ./data/\\n${g}.genome : ${g}\\n" > snpeff_db/snpEff.config
    cd snpeff_db
    snpEff -Xmx${(task.memory.toMega() * 0.75) as int}m build -gtf22 -noCheckCds -noCheckProtein -c snpEff.config ${g} > build.log 2>&1
    rm data/${g}/sequences.fa
    """
}
