/*
    Stage 5: evidence integration - resolve ARG calls into the four-tier
    confidence framework (see CLAUDE.md and issue #23).

    Built incrementally, one signal per sub-issue:
    - #24 classifier-agreement signal: DONE (CLASSIFIER_AGREEMENT below)
    - #25 replicon/mobility marker evidence: DONE (MOBILITY_MARKERS below)
    - #26 circularisation/coverage-ratio evidence: DONE (CIRCULARITY_COVERAGE below)
    - #27 validated contig-length floor: DONE (CONTIG_LENGTH_FLOOR below)
    - #28 SCCmec-cassette override: DONE (SCCMEC_OVERRIDE below)
    - #29 tier-resolution rule engine (combines all of the above): DONE (TIER_RESOLUTION below)

    Stage 5 is now complete: `tier` is the final per-ARG four-tier call.

    MOB-suite/Platon/RFPlasmid outputs stay one-file-per-sample all the way
    through subworkflows/local/plasmid_classification, so they join cleanly on
    `meta` here - no collision risk. ARG's hamronization report
    (subworkflows/local/funcscan_arg) is the one exception: it's merged across
    every sample in a single file, so it's broadcast (.first()) to every
    sample's CLASSIFIER_AGREEMENT call, which filters to its own rows itself
    (see bin/classifier_agreement.py).

    platon_tsv is joined with `remainder: true` rather than checked against
    params.plasmid_skip_platon directly - this subworkflow should only care
    whether it actually received Platon data on that channel, not re-derive
    the same decision PLASMID_CLASSIFICATION already made. Keeps this testable
    with hand-built channels regardless of global param state (an earlier
    version checked the param here too and silently ignored fixture data
    supplied straight to this subworkflow's own nf-test, since the `test`
    profile always sets plasmid_skip_platon = true).
*/

include { CLASSIFIER_AGREEMENT  } from '../../../modules/local/classifier_agreement/main'
include { MOBILITY_MARKERS      } from '../../../modules/local/mobility_markers/main'
include { CIRCULARITY_COVERAGE  } from '../../../modules/local/circularity_coverage/main'
include { CONTIG_LENGTH_FLOOR   } from '../../../modules/local/contig_length_floor/main'
include { SCCMEC_OVERRIDE       } from '../../../modules/local/sccmec_override/main'
include { TIER_RESOLUTION       } from '../../../modules/local/tier_resolution/main'

workflow EVIDENCE_INTEGRATION {
    take:
    arg_report              // path: hamronization_combined_report.tsv (all samples, stage 3)
    assembly                // channel: [ val(meta), path(fasta) ], possibly gzipped (stage 2, BACASS.out.assembly)
    mobsuite_contig_report  // channel: [ val(meta), path(contig_report.txt) ] (stage 4)
    platon_tsv              // channel: [ val(meta), path(*.tsv) ] - empty channel if --plasmid_skip_platon
    rfplasmid_prediction    // channel: [ val(meta), path(prediction.csv) ] (stage 4)
    sccmecextractor_report  // channel: [ val(meta), path(sccmec_unified_report.tsv) ] (stage 6)

    main:
    ch_versions = channel.empty()

    ch_classifier_inputs = mobsuite_contig_report
        .join(rfplasmid_prediction)
        .join(platon_tsv, remainder: true)
        .map { meta, mob, rfp, plat -> [meta, mob, rfp, plat ?: []] }

    CLASSIFIER_AGREEMENT ( ch_classifier_inputs, arg_report.first() )
    ch_versions = ch_versions.mix(CLASSIFIER_AGREEMENT.out.versions)

    MOBILITY_MARKERS ( mobsuite_contig_report, arg_report.first() )
    ch_versions = ch_versions.mix(MOBILITY_MARKERS.out.versions)

    CIRCULARITY_COVERAGE ( assembly, arg_report.first() )
    ch_versions = ch_versions.mix(CIRCULARITY_COVERAGE.out.versions)

    CONTIG_LENGTH_FLOOR ( assembly, arg_report.first() )
    ch_versions = ch_versions.mix(CONTIG_LENGTH_FLOOR.out.versions)

    SCCMEC_OVERRIDE ( sccmecextractor_report, arg_report.first() )
    ch_versions = ch_versions.mix(SCCMEC_OVERRIDE.out.versions)

    ch_tier_inputs = CLASSIFIER_AGREEMENT.out.tsv
        .join(MOBILITY_MARKERS.out.tsv)
        .join(CIRCULARITY_COVERAGE.out.tsv)
        .join(CONTIG_LENGTH_FLOOR.out.tsv)
        .join(SCCMEC_OVERRIDE.out.tsv)

    TIER_RESOLUTION ( ch_tier_inputs )
    ch_versions = ch_versions.mix(TIER_RESOLUTION.out.versions)

    emit:
    classifier_agreement = CLASSIFIER_AGREEMENT.out.tsv  // channel: [ val(meta), path(*.classifier_agreement.tsv) ]
    mobility_markers      = MOBILITY_MARKERS.out.tsv     // channel: [ val(meta), path(*.mobility_markers.tsv) ]
    circularity_coverage  = CIRCULARITY_COVERAGE.out.tsv // channel: [ val(meta), path(*.circularity_coverage.tsv) ]
    contig_length_floor   = CONTIG_LENGTH_FLOOR.out.tsv  // channel: [ val(meta), path(*.contig_length_floor.tsv) ]
    sccmec_override       = SCCMEC_OVERRIDE.out.tsv      // channel: [ val(meta), path(*.sccmec_override.tsv) ]
    tier_resolution        = TIER_RESOLUTION.out.tsv      // channel: [ val(meta), path(*.tier_resolution.tsv) ] - the final per-ARG four-tier call (stage 7, #22, validates this)
    versions              = ch_versions
}
