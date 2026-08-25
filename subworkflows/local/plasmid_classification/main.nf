/*
    Run plasmid classification tools (stage 4): MOB-suite, Platon, RFPlasmid.

    CLAUDE.md frames stage 4 as three classifiers run in parallel, whose
    agreement/disagreement stage 5 (evidence integration, not yet built) reads
    directly - this subworkflow is where all three eventually converge.

    All three classifiers (MOB-suite #8, Platon #9, RFPlasmid #10) are wired
    in - stage 4 is complete.

    mob_recon's and RFPLASMID's own scripts decompress a gzipped assembly
    internally, and PLATON (modules/local/platon) copies that same handling,
    unlike spaTyper/staphopia-sccmec in subworkflows/local/typing - no
    separate GUNZIP step is needed here for any of the three.

    Platon and RFPlasmid both need an external database download (2.8 GB and
    ~645 MB respectively) before they can run at all - too heavy for routine
    dev-loop or CI runs, matching how subworkflows/local/funcscan_arg/main.nf
    lets RGI/DeepARG be skipped. Both default to running in real pipeline
    runs and are only skipped in this subworkflow's own nf-test.
*/

include { MOBSUITE_RECON      } from '../../../modules/nf-core/mobsuite/recon/main'
include { PLATON              } from '../../../modules/local/platon/main'
include { RFPLASMID           } from '../../../modules/local/rfplasmid/main'
include { RFPLASMID_DOWNLOADDB } from '../../../modules/local/rfplasmid/downloaddb/main'
include { UNTAR as UNTAR_PLATON_DB } from '../../../modules/nf-core/untar/main'

workflow PLASMID_CLASSIFICATION {
    take:
    assembly // tuple val(meta), path(fasta), possibly gzipped

    main:
    ch_versions = channel.empty()

    MOBSUITE_RECON ( assembly )

    // Platon needs a ~2.8 GB reference database (Zenodo download, see below).
    // plasmid_skip_platon defaults to false (it runs in real pipeline runs).
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

    // RFPlasmid needs two DIAMOND databases built from ~645 MB of protein
    // FASTA (modules/local/rfplasmid/downloaddb) - the Random Forest models
    // themselves ship bundled in the container, but these do not.
    // plasmid_skip_rfplasmid defaults to false (it runs in real pipeline runs).
    ch_rfplasmid_prediction      = channel.empty()
    ch_rfplasmid_prediction_full = channel.empty()
    if (!params.plasmid_skip_rfplasmid) {
        if (params.plasmid_rfplasmid_db) {
            ch_rfplasmid_db = channel.fromPath("${params.plasmid_rfplasmid_db}/*.dmnd", checkIfExists: true).collect()
        }
        else {
            RFPLASMID_DOWNLOADDB ()
            ch_versions = ch_versions.mix(RFPLASMID_DOWNLOADDB.out.versions)
            ch_rfplasmid_db = RFPLASMID_DOWNLOADDB.out.db.collect()
        }

        RFPLASMID ( assembly, params.plasmid_rfplasmid_species, ch_rfplasmid_db )
        ch_rfplasmid_prediction      = RFPLASMID.out.prediction
        ch_rfplasmid_prediction_full = RFPLASMID.out.prediction_full
    }

    emit:
    mobsuite_contig_report    = MOBSUITE_RECON.out.contig_report    // channel: [ val(meta), path(contig_report.txt) ] - per-contig chromosome/plasmid-group call
    mobsuite_chromosome       = MOBSUITE_RECON.out.chromosome       // channel: [ val(meta), path(chromosome.fasta) ]
    mobsuite_plasmids         = MOBSUITE_RECON.out.plasmids         // channel: [ val(meta), path(plasmid_*.fasta) ]
    mobsuite_mobtyper_results = MOBSUITE_RECON.out.mobtyper_results // channel: [ val(meta), path(mobtyper_results.txt) ] - replicon/relaxase typing per plasmid, feeds stage 5b (#25)
    platon_tsv                = ch_platon_tsv                      // channel: [ val(meta), path(*.tsv) ] - per-contig chromosome/plasmid call
    platon_json                = ch_platon_json                    // channel: [ val(meta), path(*.json) ]
    rfplasmid_prediction       = ch_rfplasmid_prediction           // channel: [ val(meta), path(prediction.csv) ] - per-contig chromosome/plasmid call
    rfplasmid_prediction_full  = ch_rfplasmid_prediction_full      // channel: [ val(meta), path(prediction_full.csv) ]
    versions                   = ch_versions
}
