// Not in the nf-core/modules registry - wraps bin/tier_resolution.py
// (stage 5f, issue #29). Same mulled python+pandas container as the other
// evidence_integration modules. Unlike them, this one only needs the five
// already-per-sample signal tsvs (each already filtered to one sample by
// the upstream modules) - no hamronization/sample-id broadcast needed here.
process TIER_RESOLUTION {
    tag "$meta.id"
    label 'process_low'

    conda "conda-forge::python=3.11.0 conda-forge::biopython=1.80 conda-forge::pandas=1.5.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' :
        'biocontainers/mulled-v2-27978155697a3671f3ef9aead4b5c823a02cc0b7:548df772fe13c0232a7eab1bc1deb98b495a05ab-0' }"

    input:
    tuple val(meta), path(classifier_agreement), path(mobility_markers), path(circularity_coverage), path(contig_length_floor), path(sccmec_override)

    output:
    tuple val(meta), path("*.tier_resolution.tsv"), emit: tsv
    path "versions.yml"                            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    tier_resolution.py \\
        --classifier-agreement ${classifier_agreement} \\
        --mobility-markers ${mobility_markers} \\
        --circularity-coverage ${circularity_coverage} \\
        --contig-length-floor ${contig_length_floor} \\
        --sccmec-override ${sccmec_override} \\
        --output ${prefix}.tier_resolution.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        tier_resolution: \$(tier_resolution.py --version | sed 's/tier_resolution //g')
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
