// Not in the nf-core/modules registry - wraps bin/circularity_coverage.py
// (stage 5c, issue #26). Same mulled python+pandas container and
// sample-broadcast/filter pattern as modules/local/classifier_agreement -
// see that module's and bin/classifier_agreement.py's comments for why.
process CIRCULARITY_COVERAGE {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), path(assembly)
    path(hamronization_tsv)

    output:
    tuple val(meta), path("*.circularity_coverage.tsv"), emit: tsv
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    circularity_coverage.py \\
        --hamronization ${hamronization_tsv} \\
        --sample-id ${meta.id} \\
        --assembly ${assembly} \\
        --coverage-ratio-threshold ${params.plasmid_coverage_ratio_threshold} \\
        --output ${prefix}.circularity_coverage.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        circularity_coverage: \$(circularity_coverage.py --version | sed 's/circularity_coverage //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
