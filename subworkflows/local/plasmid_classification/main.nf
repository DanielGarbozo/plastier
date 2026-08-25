/*
    Run plasmid classification tools (stage 4): MOB-suite, Platon, RFPlasmid.

    CLAUDE.md frames stage 4 as three classifiers run in parallel, whose
    agreement/disagreement stage 5 (evidence integration, not yet built) reads
    directly - this subworkflow is where all three eventually converge.

    MOB-suite (issue #8) and Platon (issue #9) are wired in. RFPlasmid
    (issue #10) is not yet vendored - add it the same way: include the
    module, call it on `assembly`, mix its versions in, and emit its
    per-contig classification output alongside the other two.

    mob_recon's own script decompresses a gzipped assembly internally, and
    PLATON (modules/local/platon) copies that same handling, unlike
    spaTyper/staphopia-sccmec in subworkflows/local/typing - no separate
    GUNZIP step is needed here for either.
*/

include { MOBSUITE_RECON } from '../../../modules/nf-core/mobsuite/recon/main'
include { PLATON         } from '../../../modules/local/platon/main'
include { UNTAR as UNTAR_PLATON_DB } from '../../../modules/nf-core/untar/main'

workflow PLASMID_CLASSIFICATION {
    take:
    assembly // tuple val(meta), path(fasta), possibly gzipped

    main:
    ch_versions = channel.empty()

    MOBSUITE_RECON ( assembly )

    // Platon needs a ~2.8 GB reference database (Zenodo download, see below) -
    // matching how subworkflows/local/funcscan_arg/main.nf lets RGI/DeepARG be
    // skipped (arg_skip_rgi etc.) since their databases are too heavy for
    // routine dev-loop or CI runs. plasmid_skip_platon defaults to false (it
    // runs in real pipeline runs) and is only set true in this subworkflow's
    // own nf-test.
    ch_platon_tsv  = channel.empty()
    ch_platon_json = channel.empty()
    if (!params.plasmid_skip_platon) {
        if (params.plasmid_platon_db) {
            ch_platon_db = file(params.plasmid_platon_db, checkIfExists: true)
        }
        else {
            UNTAR_PLATON_DB ( [ [], file('https://zenodo.org/record/4066768/files/db.tar.gz', checkIfExists: true) ] )
            ch_platon_db = UNTAR_PLATON_DB.out.untar.map { it[1] }
        }

        PLATON ( assembly, ch_platon_db )
        ch_platon_tsv  = PLATON.out.tsv
        ch_platon_json = PLATON.out.json
    }

    emit:
    mobsuite_contig_report    = MOBSUITE_RECON.out.contig_report    // channel: [ val(meta), path(contig_report.txt) ] - per-contig chromosome/plasmid-group call
    mobsuite_chromosome       = MOBSUITE_RECON.out.chromosome       // channel: [ val(meta), path(chromosome.fasta) ]
    mobsuite_plasmids         = MOBSUITE_RECON.out.plasmids         // channel: [ val(meta), path(plasmid_*.fasta) ]
    mobsuite_mobtyper_results = MOBSUITE_RECON.out.mobtyper_results // channel: [ val(meta), path(mobtyper_results.txt) ] - replicon/relaxase typing per plasmid, feeds stage 5b (#25)
    platon_tsv                = ch_platon_tsv                       // channel: [ val(meta), path(*.tsv) ] - per-contig chromosome/plasmid call
    platon_json                = ch_platon_json                     // channel: [ val(meta), path(*.json) ]
    versions                  = ch_versions
}
