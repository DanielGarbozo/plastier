/*
    Run ARG screening tools

    Vendored from nf-core/funcscan 4.0.0 (subworkflows/local/arg.nf).
    Only the ARG (antimicrobial resistance gene) screening path is vendored here -
    AMP/BGC/CAZyme screening and taxonomic classification are out of scope for plastier
    and were never brought in.

    NOTE: deviation from upstream funcscan 4.0.0 - upstream's ARG workflow only
    `emit: versions`, leaving the merged hAMRonization report to be written to disk
    solely via HAMRONIZATION_SUMMARIZE's own publishDir. plastier's stage 5 (evidence
    integration, not yet built) needs to programmatically consume that merged report,
    so a second output - `emit: report` - has been added below, threading through
    HAMRONIZATION_SUMMARIZE.out.tsv. HAMRONIZATION_SUMMARIZE is called from inside
    this file both upstream and here (it is not called from the parent funcscan.nf),
    so no call needed to be moved in from the parent workflow.
*/

include { ABRICATE_RUN                     } from '../../../modules/nf-core/abricate/run'
include { AMRFINDERPLUS_UPDATE             } from '../../../modules/nf-core/amrfinderplus/update'
include { AMRFINDERPLUS_RUN                } from '../../../modules/nf-core/amrfinderplus/run'
include { DEEPARG_DOWNLOADDATA             } from '../../../modules/nf-core/deeparg/downloaddata'
include { DEEPARG_PREDICT                  } from '../../../modules/nf-core/deeparg/predict'
include { FARGENE                          } from '../../../modules/nf-core/fargene'
include { GUNZIP as ARG_FARGENE_GUNZIP     } from '../../../modules/nf-core/gunzip/main'
include { RGI_CARDANNOTATION               } from '../../../modules/nf-core/rgi/cardannotation'
include { RGI_MAIN                         } from '../../../modules/nf-core/rgi/main'
include { UNTAR as UNTAR_CARD              } from '../../../modules/nf-core/untar'
include { TABIX_BGZIP as ARG_TABIX_BGZIP   } from '../../../modules/nf-core/tabix/bgzip'
include { MERGE_TAXONOMY_HAMRONIZATION     } from '../../../modules/local/funcscan_arg/merge_taxonomy_hamronization'
include { HAMRONIZATION_RGI                } from '../../../modules/nf-core/hamronization/rgi'
include { HAMRONIZATION_FARGENE            } from '../../../modules/nf-core/hamronization/fargene'
include { HAMRONIZATION_SUMMARIZE          } from '../../../modules/nf-core/hamronization/summarize'
include { HAMRONIZATION_ABRICATE           } from '../../../modules/nf-core/hamronization/abricate'
include { HAMRONIZATION_DEEPARG            } from '../../../modules/nf-core/hamronization/deeparg'
include { HAMRONIZATION_AMRFINDERPLUS      } from '../../../modules/nf-core/hamronization/amrfinderplus'
include { ARGNORM as ARGNORM_DEEPARG       } from '../../../modules/nf-core/argnorm'
include { ARGNORM as ARGNORM_ABRICATE      } from '../../../modules/nf-core/argnorm'
include { ARGNORM as ARGNORM_AMRFINDERPLUS } from '../../../modules/nf-core/argnorm'

workflow ARG {
    take:
    fastas // tuple val(meta), path(contigs)
    annotations
    tsvs // tuple val(meta), path(MMSEQS_CREATETSV.out.tsv)

    main:
    ch_versions = channel.empty()

    // Prepare HAMRONIZATION reporting channel
    ch_input_to_hamronization_summarize = channel.empty()

    // AMRfinderplus run
    // Prepare channel for database
    if (!params.arg_skip_amrfinderplus && params.arg_amrfinderplus_db) {
        ch_amrfinderplus_db = channel.fromPath(params.arg_amrfinderplus_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_amrfinderplus && !params.arg_amrfinderplus_db) {
        AMRFINDERPLUS_UPDATE()
        // NOTE: deviation from upstream - the vendored modules/nf-core/amrfinderplus/update
        // emits software versions via the `versions` topic instead of a classic
        // `.out.versions` channel, so there is nothing to mix in here anymore.
        ch_amrfinderplus_db = AMRFINDERPLUS_UPDATE.out.db
    }

    if (!params.arg_skip_amrfinderplus) {
        AMRFINDERPLUS_RUN(fastas, ch_amrfinderplus_db)
        // NOTE: deviation from upstream - modules/nf-core/amrfinderplus/run emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.

        // Reporting
        HAMRONIZATION_AMRFINDERPLUS(AMRFINDERPLUS_RUN.out.report, 'tsv', AMRFINDERPLUS_RUN.out.tool_version, AMRFINDERPLUS_RUN.out.db_version)
        // NOTE: deviation from upstream - modules/nf-core/hamronization/amrfinderplus emits
        // software versions via the `versions` topic instead of a classic `.out.versions`
        // channel.
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_AMRFINDERPLUS.out.tsv)

        if (!params.arg_skip_argnorm) {
            ch_input_to_argnorm_amrfinderplus = HAMRONIZATION_AMRFINDERPLUS.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_AMRFINDERPLUS(ch_input_to_argnorm_amrfinderplus, 'amrfinderplus', 'ncbi')
            // NOTE: deviation from upstream - modules/nf-core/argnorm emits software versions
            // via the `versions` topic instead of a classic `.out.versions` channel.
        }
    }

    // fARGene run
    if (!params.arg_skip_fargene) {
        ch_fargene_classes = channel.fromList(params.arg_fargene_hmmmodel.tokenize(','))

        // Unlike AMRfinderplus/ABRicate/RGI, fARGene does not auto-decompress gzip input
        // and exits with "input file(s) must be FASTA" when given one - same problem
        // TYPING_GUNZIP already solves for the typing subworkflow, mirrored here.
        fastas
            .branch { _meta, fasta ->
                gzip:     "${fasta}".endsWith('.gz')
                not_gzip: true
            }
            .set { ch_fastas_for_fargene_gunzip }

        ARG_FARGENE_GUNZIP ( ch_fastas_for_fargene_gunzip.gzip )

        ch_fastas_uncompressed = ARG_FARGENE_GUNZIP.out.gunzip.mix(ch_fastas_for_fargene_gunzip.not_gzip)

        ch_fargene_input = ch_fastas_uncompressed
            .combine(ch_fargene_classes)
            .map { meta, fastafiles, hmm_class ->
                def meta_new = meta.clone()
                meta_new['hmm_class'] = hmm_class
                [meta_new, fastafiles, hmm_class]
            }
            .multiMap {
                fastas: [it[0], it[1]]
                hmmclass: it[2]
            }

        FARGENE(ch_fargene_input.fastas, ch_fargene_input.hmmclass)
        // NOTE: deviation from upstream - modules/nf-core/fargene emits software versions via
        // the `versions` topic instead of a classic `.out.versions` channel.

        // Reporting
        // Note: currently hardcoding versions, has to be updated with every fARGene-update
        HAMRONIZATION_FARGENE(FARGENE.out.hmm_genes.transpose(), 'tsv', '0.1', '0.1')
        // NOTE: deviation from upstream - modules/nf-core/hamronization/fargene emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_FARGENE.out.tsv)
    }

    // RGI run
    if (!params.arg_skip_rgi) {

        if (!params.arg_rgi_db) {

            // Download and untar CARD
            UNTAR_CARD([[], file('assets/card.tar.bz2', checkIfExists: true)])
            // NOTE: deviation from upstream - modules/nf-core/untar emits software versions via
            // the `versions` topic instead of a classic `.out.versions` channel.
            rgi_db = UNTAR_CARD.out.untar.map { it[1] }
            RGI_CARDANNOTATION(rgi_db)
            card = RGI_CARDANNOTATION.out.db
            ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
        }
        else {

            // Use user-supplied database
            rgi_db = file(params.arg_rgi_db, checkIfExists: true)
            if (!rgi_db.contains("card_database_processed")) {
                RGI_CARDANNOTATION(rgi_db)
                card = RGI_CARDANNOTATION.out.db
                ch_versions = ch_versions.mix(RGI_CARDANNOTATION.out.versions)
            }
            else {
                card = rgi_db
            }
        }

        RGI_MAIN(fastas, card, [])
        ch_versions = ch_versions.mix(RGI_MAIN.out.versions)

        // Reporting
        HAMRONIZATION_RGI(RGI_MAIN.out.tsv, 'tsv', RGI_MAIN.out.tool_version, RGI_MAIN.out.db_version)
        // NOTE: deviation from upstream - modules/nf-core/hamronization/rgi emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_RGI.out.tsv)
    }

    // DeepARG prepare download
    if (!params.arg_skip_deeparg && params.arg_deeparg_db) {
        ch_deeparg_db = channel.fromPath(params.arg_deeparg_db, checkIfExists: true)
            .first()
    }
    else if (!params.arg_skip_deeparg && !params.arg_deeparg_db) {
        DEEPARG_DOWNLOADDATA()
        // NOTE: deviation from upstream - modules/nf-core/deeparg/downloaddata emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.
        ch_deeparg_db = DEEPARG_DOWNLOADDATA.out.db
    }

    // DeepARG run
    if (!params.arg_skip_deeparg) {

        annotations
            .map { it ->
                def meta = it[0]
                def anno = it[1]
                def model = params.arg_deeparg_model

                [meta, anno, model]
            }
            .set { ch_input_for_deeparg }

        DEEPARG_PREDICT(ch_input_for_deeparg, ch_deeparg_db)
        // NOTE: deviation from upstream - modules/nf-core/deeparg/predict emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.

        // Reporting
        // Note: currently hardcoding versions as unreported by DeepARG
        // Make sure to update on version bump.
        ch_input_to_hamronization_deeparg = DEEPARG_PREDICT.out.arg.mix(DEEPARG_PREDICT.out.potential_arg)
        HAMRONIZATION_DEEPARG(ch_input_to_hamronization_deeparg, 'tsv', '1.0.4', params.arg_deeparg_db_version)
        // NOTE: deviation from upstream - modules/nf-core/hamronization/deeparg emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_DEEPARG.out.tsv)

        if (!params.arg_skip_argnorm) {
            ch_input_to_argnorm_deeparg = HAMRONIZATION_DEEPARG.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_DEEPARG(ch_input_to_argnorm_deeparg, 'deeparg', 'deeparg')
            // NOTE: deviation from upstream - modules/nf-core/argnorm emits software versions
            // via the `versions` topic instead of a classic `.out.versions` channel.
        }
    }

    // ABRicate run
    if (!params.arg_skip_abricate) {
        abricate_dbdir = params.arg_abricate_db ? file(params.arg_abricate_db, checkIfExists: true) : []
        ABRICATE_RUN(fastas, abricate_dbdir)
        // NOTE: deviation from upstream - modules/nf-core/abricate/run emits software versions
        // via the `versions` topic instead of a classic `.out.versions` channel.

        HAMRONIZATION_ABRICATE(ABRICATE_RUN.out.report, 'tsv', '1.0.1', '2021-Mar-27')
        // NOTE: deviation from upstream - modules/nf-core/hamronization/abricate emits software
        // versions via the `versions` topic instead of a classic `.out.versions` channel.
        ch_input_to_hamronization_summarize = ch_input_to_hamronization_summarize.mix(HAMRONIZATION_ABRICATE.out.tsv)

        if ((params.arg_abricate_db_id == 'ncbi' || params.arg_abricate_db_id == 'resfinder' || params.arg_abricate_db_id == 'argannot' || params.arg_abricate_db_id == 'megares') && !params.arg_skip_argnorm) {
            ch_input_to_argnorm_abricate = HAMRONIZATION_ABRICATE.out.tsv.filter { meta, file -> !file.isEmpty() }
            ARGNORM_ABRICATE(ch_input_to_argnorm_abricate, 'abricate', params.arg_abricate_db_id)
            // NOTE: deviation from upstream - modules/nf-core/argnorm emits software versions
            // via the `versions` topic instead of a classic `.out.versions` channel.
        }
    }

    ch_input_to_hamronization_summarize
        .map {
            it[1]
        }
        .collect()
        .set { ch_input_for_hamronization_summarize }

    HAMRONIZATION_SUMMARIZE(ch_input_for_hamronization_summarize, params.arg_hamronization_summarizeformat)
    // NOTE: deviation from upstream - modules/nf-core/hamronization/summarize emits software
    // versions via the `versions` topic instead of a classic `.out.versions` channel.

    // MERGE_TAXONOMY
    if (params.run_taxa_classification) {

        ch_mmseqs_taxonomy_list = tsvs.map { it[1] }.collect()
        MERGE_TAXONOMY_HAMRONIZATION(HAMRONIZATION_SUMMARIZE.out.tsv, ch_mmseqs_taxonomy_list)
        ch_versions = ch_versions.mix(MERGE_TAXONOMY_HAMRONIZATION.out.versions)

        ch_tabix_input = channel.of(['id': 'hamronization_combined_report'])
            .combine(MERGE_TAXONOMY_HAMRONIZATION.out.tsv)

        ARG_TABIX_BGZIP(ch_tabix_input)
        // NOTE: deviation from upstream - modules/nf-core/tabix/bgzip emits software versions
        // via the `versions` topic instead of a classic `.out.versions` channel.
    }

    emit:
    versions = ch_versions
    // NOTE: deviation from upstream funcscan 4.0.0 - upstream ARG only `emit: versions`.
    // plastier's stage 5 (evidence integration, future work) needs to programmatically
    // consume the merged hAMRonization report rather than relying on its publishDir side
    // effect, so it is threaded through here as a second output.
    report = HAMRONIZATION_SUMMARIZE.out.tsv // path(hamronization_combined_report.tsv)
}
