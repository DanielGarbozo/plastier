// Not in the nf-core/modules registry - hand-written to match nf-core module
// conventions (see modules/local/platon/main.nf, the closest local precedent).
//
// Confirmed by inspecting the actual biocontainers image (0.0.18--pyhdfd78af_0):
// no --db/--initialize flag exists in this version, and each genus's trained
// Random Forest model ships bundled in the package (e.g. Staphylococcus.rfo).
// But two DIAMOND databases (plasmiddb_cge.dmnd, plasmiddb_total.dmnd) are
// NOT bundled and are required - rfplasmid looks for them at a fixed path
// next to its own installation (no flag to redirect them), so this module
// resolves that path at runtime and copies in whatever
// modules/local/rfplasmid/downloaddb produced.
//
// rfplasmid.py globs `*.fasta` inside --input, which must be a directory, and
// os.chdir()s into --out before writing prediction.csv/prediction_full.csv -
// see classification.R in the package for the exact column layout
// (prediction, votes chromosomal, votes plasmid, + feature columns).
process RFPLASMID {
    tag "$meta.id"
    label 'process_medium'
    // rfplasmid has no flag to redirect its .dmnd lookup path (see comment
    // above) - the *.dmnd files have to be copied into its own install
    // directory, which the default host-user-mapped container (docker
    // .runOptions '-u $(id -u):$(id -g)' in nextflow.config) cannot write to
    // (that directory is root-owned inside the image). Root only for this
    // process, purely to satisfy that write.
    containerOptions "-u root"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rfplasmid:0.0.18--pyhdfd78af_0' :
        'quay.io/biocontainers/rfplasmid:0.0.18--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)
    val species
    path db // the two *.dmnd files from RFPLASMID_DOWNLOADDB (or a user-supplied equivalent)

    output:
    tuple val(meta), path("results/prediction.csv")      , emit: prediction // per-contig chromosome/plasmid call - feeds stage 5a (#24)
    tuple val(meta), path("results/prediction_full.csv") , emit: prediction_full, optional: true
    tuple val("${task.process}"), val('rfplasmid'), eval('rfplasmid --version 2>&1 | sed "s/^.*version //"'), emit: versions_rfplasmid, topic: versions

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

    # rfplasmid.py globs `*.fasta` inside --input, so the staged copy must use
    # that exact extension regardless of the original file's name.
    mkdir -p input_dir
    cp ${fasta_name} input_dir/${prefix}.fasta

    RFPLASMID_PKG_DIR=\$(python3 -c "import RFPlasmid, os; print(os.path.dirname(RFPlasmid.__file__))")
    cp *.dmnd \${RFPLASMID_PKG_DIR}/

    rfplasmid \\
        --species ${species} \\
        --input input_dir \\
        --out results \\
        --threads ${task.cpus} \\
        --jelly \\
        ${args}
    """

    stub:
    """
    mkdir -p results
    touch results/prediction.csv
    touch results/prediction_full.csv
    """
}
