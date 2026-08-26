// Not in the nf-core/modules registry (no bioconda package exists yet) -
// hand-written to match nf-core module conventions. pip-installable
// (PyPI: sccmecextractor) with a pinned Docker image published by the tool's
// own author (alisonmacfadyen/sccmecextractor) - verified a versioned tag
// exists (v1.6.0, matching the current PyPI release) rather than depending
// on :latest.
//
// Chosen over staphopia-sccmec (already in subworkflows/local/typing) for
// stage 5e (issue #28) specifically because staphopia-sccmec only reports
// whether a genome has each SCCmec type present/absent (confirmed by
// inspecting its real --json output, which carries zero positional
// information) - it cannot tell you WHERE the cassette is, which stage 5e's
// override needs to test ARG-coordinate overlap. SCCmecExtractor locates the
// actual att-site boundaries per contig.
process SCCMECEXTRACTOR {
    tag "$meta.id"
    label 'process_medium'

    // No bioconda package exists for this tool, so unlike every other module
    // in this repo there is no galaxyproject-hosted singularity mirror to
    // point at - singularity/apptainer pull the same Docker Hub image
    // directly instead (`docker://` prefix), same tag either way. Both
    // branches spell out `docker.io/` explicitly - this repo's nextflow.config
    // sets docker.registry = 'quay.io' as the default, which silently
    // prepends to any unqualified image name and 401s trying to resolve this
    // image there (found by running the actual test, not by inspection).
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://docker.io/alisonmacfadyen/sccmecextractor:v1.6.0' :
        'docker.io/alisonmacfadyen/sccmecextractor:v1.6.0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("results/sccmec_unified_report.tsv"), emit: report
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def is_compressed = fasta.getName().endsWith(".gz") ? true : false
    def fasta_name = fasta.getName().replace(".gz", "")
    """
    if [ "${is_compressed}" == "true" ]; then
        gzip -c -d ${fasta} > ${fasta_name}
    fi

    sccmec-pipeline \\
        -f ${fasta_name} \\
        -o results \\
        -t ${task.cpus} \\
        ${args}

    # sccmec-pipeline has no --version flag (it prints its usage/help text
    # instead, which broke versions.yml's YAML when captured directly) -
    # hardcoded to match the pinned container tag, same as done for fARGene.
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sccmecextractor: "1.6.0"
    END_VERSIONS
    """
}
