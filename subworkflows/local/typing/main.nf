/*
    Run strain typing tools (stage 6): MLST, spa typing, SCCmec typing.

    All three run directly off BACASS's assembly - no re-annotation needed.
    spaTyper and staphopia-sccmec are S. aureus-specific; MLST auto-detects
    the matching PubMLST scheme from the assembly itself.

    BACASS's default short-read assembler (unicycler) emits a gzipped
    scaffolds fasta. mlst reads gzipped fastas natively (its own upstream
    module test feeds it one directly), but spaTyper and staphopia-sccmec do
    not - they try to decode the gzip magic bytes as text and crash with a
    UnicodeDecodeError. Only those two get a decompressed copy.

    Typing output is not yet consumed anywhere - stage 5 (evidence integration,
    not yet built) is what will use a typed SCCmec cassette's coordinates to
    resolve an ARG sitting inside it to the chromosomal tier instead of ambiguous,
    per the four-tier framework in CLAUDE.md.
*/

include { GUNZIP as TYPING_GUNZIP } from '../../../modules/nf-core/gunzip/main'
include { MLST                    } from '../../../modules/nf-core/mlst/main'
include { SPATYPER                } from '../../../modules/nf-core/spatyper/main'
include { STAPHOPIASCCMEC         } from '../../../modules/nf-core/staphopiasccmec/main'

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

    emit:
    mlst_tsv   = MLST.out.tsv            // channel: [ val(meta), path(tsv) ]
    spa_tsv    = SPATYPER.out.tsv        // channel: [ val(meta), path(tsv) ]
    sccmec_tsv = STAPHOPIASCCMEC.out.tsv // channel: [ val(meta), path(tsv) ]
    versions   = ch_versions
}
