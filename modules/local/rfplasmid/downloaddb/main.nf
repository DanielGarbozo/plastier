// Not in the nf-core/modules registry - hand-written to match nf-core module
// conventions. Mirrors modules/local/rfplasmid/main.nf's getdb.sh, which the
// upstream package ships but never invokes automatically: two protein FASTAs
// from https://klif.uu.nl/download/plasmid_db/ (a small ~300 KB "cge" set and
// a ~645 MB "total" set) get indexed into DIAMOND databases. rfplasmid looks
// for both *.dmnd files at a fixed path next to its own installation, not a
// user-configurable one - see modules/local/rfplasmid/main.nf for how the
// output of this process gets copied into place before each prediction.
process RFPLASMID_DOWNLOADDB {
    label 'process_low'

    conda "${moduleDir}/../environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rfplasmid:0.0.18--pyhdfd78af_0' :
        'quay.io/biocontainers/rfplasmid:0.0.18--pyhdfd78af_0' }"

    output:
    path("*.dmnd")     , emit: db
    path("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget -O plasmiddb_cge.faa https://klif.uu.nl/download/plasmid_db/plasmiddb_cge.faa
    wget -O plasmiddb_total.faa https://klif.uu.nl/download/plasmid_db/plasmiddb_total.faa

    diamond makedb --in plasmiddb_cge.faa -d plasmiddb_cge
    diamond makedb --in plasmiddb_total.faa -d plasmiddb_total

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond version 2>&1 | sed 's/^diamond version //')
    END_VERSIONS
    """
}
