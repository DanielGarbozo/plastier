/*
    Run strain typing tools (stage 6): MLST, spa typing, SCCmec typing
    (presence/absence) + SCCmec cassette extraction (boundary coordinates).

    All four run directly off BACASS's assembly - no re-annotation needed.
    spaTyper, staphopia-sccmec, and SCCmecExtractor are S. aureus-specific;
    MLST auto-detects the matching PubMLST scheme from the assembly itself.

    BACASS's default short-read assembler (unicycler) emits a gzipped
    scaffolds fasta. mlst reads gzipped fastas natively (its own upstream
    module test feeds it one directly), but spaTyper and staphopia-sccmec do
    not - they try to decode the gzip magic bytes as text and crash with a
    UnicodeDecodeError. Only those two get a decompressed copy.
    SCCMECEXTRACTOR handles its own decompression internally (see
    modules/local/sccmecextractor), matching the pattern used by
    subworkflows/local/plasmid_classification's tools.

    Two different SCCmec tools, on purpose: staphopia-sccmec only reports
    whether the genome has each SCCmec type present/absent (confirmed against
    its real --json output - no positional information at all), which is
    enough for typing but not for stage 5e (issue #28)'s override, which
    needs to test whether an ARG's coordinates actually fall inside the
    cassette. SCCmecExtractor locates the real att-site boundaries per
    contig - see subworkflows/local/evidence_integration for how that gets
    used.
*/

include { GUNZIP as TYPING_GUNZIP } from '../../../modules/nf-core/gunzip/main'
include { MLST                    } from '../../../modules/nf-core/mlst/main'
include { SPATYPER                } from '../../../modules/nf-core/spatyper/main'
include { STAPHOPIASCCMEC         } from '../../../modules/nf-core/staphopiasccmec/main'
include { SCCMECEXTRACTOR         } from '../../../modules/local/sccmecextractor/main'

workflow TYPING {
    take:
    assembly // tuple val(meta), path(fasta), possibly gzipped

    main:
    ch_versions = channel.empty()

    MLST ( assembly )

    assembly
        .branch { _meta, fasta ->
            gzip:     "${fasta}".endsWith('.gz')
            not_gzip: true
        }
        .set { ch_assembly_for_gunzip }

    TYPING_GUNZIP ( ch_assembly_for_gunzip.gzip )

    ch_assembly_uncompressed = TYPING_GUNZIP.out.gunzip.mix(ch_assembly_for_gunzip.not_gzip)

    SPATYPER (
        ch_assembly_uncompressed,
        [],
        []
    )
    ch_versions = ch_versions.mix(SPATYPER.out.versions)

    STAPHOPIASCCMEC ( ch_assembly_uncompressed )

    SCCMECEXTRACTOR ( assembly )
    ch_versions = ch_versions.mix(SCCMECEXTRACTOR.out.versions)

    emit:
    mlst_tsv             = MLST.out.tsv               // channel: [ val(meta), path(tsv) ]
    spa_tsv              = SPATYPER.out.tsv           // channel: [ val(meta), path(tsv) ]
    sccmec_tsv           = STAPHOPIASCCMEC.out.tsv    // channel: [ val(meta), path(tsv) ]
    sccmecextractor_report = SCCMECEXTRACTOR.out.report // channel: [ val(meta), path(sccmec_unified_report.tsv) ] - feeds stage 5e (#28)
    versions             = ch_versions
}
