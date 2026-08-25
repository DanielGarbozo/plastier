// Not in the nf-core/modules registry - wraps bin/classifier_agreement.py
// (stage 5a, issue #24). Reuses the same mulled python+pandas container as
// modules/local/funcscan_arg/merge_taxonomy_hamronization, already proven in
// this repo.
//
// `hamronization_tsv` is the single report merged across every sample in the
// run (subworkflows/local/funcscan_arg's HAMRONIZATION_SUMMARIZE collects all
// samples before writing it) - the calling subworkflow broadcasts the same
// path to every sample's invocation of this process, and the script filters
// to `meta.id`'s rows itself via --sample-id.
process CLASSIFIER_AGREEMENT {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), path(mobsuite_contig_report), path(rfplasmid_prediction), path(platon_tsv)
    path(hamronization_tsv)

    output:
    tuple val(meta), path("*.classifier_agreement.tsv"), emit: tsv
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def platon_arg = platon_tsv ? "--platon ${platon_tsv}" : ''
    """
    classifier_agreement.py \\
        --hamronization ${hamronization_tsv} \\
        --sample-id ${meta.id} \\
        --mobsuite ${mobsuite_contig_report} \\
        --rfplasmid ${rfplasmid_prediction} \\
        ${platon_arg} \\
        --output ${prefix}.classifier_agreement.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        classifier_agreement: \$(classifier_agreement.py --version | sed 's/classifier_agreement //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
