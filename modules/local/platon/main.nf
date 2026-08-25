// Not in the nf-core/modules registry - hand-written to match nf-core module
// conventions (see modules/nf-core/mobsuite/recon/main.nf for the closest
// registry equivalent this was modelled on).
process PLATON {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/platon:1.8--pyhdfd78af_0' :
        'quay.io/biocontainers/platon:1.8--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    path db

    output:
    tuple val(meta), path("*.tsv")             , emit: tsv // per-contig chromosome/plasmid call - feeds stage 5a (#24)
    tuple val(meta), path("*.json")            , emit: json,       optional: true
    tuple val(meta), path("*.chromosome.fasta"), emit: chromosome, optional: true
    tuple val(meta), path("*.plasmid.fasta")   , emit: plasmids,   optional: true
    tuple val("${task.process}"), val('platon'), eval('platon --version 2>&1 | sed "s/^.*platon //"'), emit: versions_platon, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    platon \\
        --db ${db} \\
        --prefix ${prefix} \\
        --output . \\
        --threads ${task.cpus} \\
        ${args} \\
        ${fasta_name}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.json
    """
}
