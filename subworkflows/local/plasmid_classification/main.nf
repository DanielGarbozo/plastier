/*
    Run plasmid classification tools (stage 4): MOB-suite, Platon, RFPlasmid.

    CLAUDE.md frames stage 4 as three classifiers run in parallel, whose
    agreement/disagreement stage 5 (evidence integration, not yet built) reads
    directly - this subworkflow is where all three eventually converge.

    Only MOB-suite (issue #8) is wired in so far. Platon (issue #9) and
    RFPlasmid (issue #10) are not yet vendored - add each the same way: include
    the module, call it on `assembly`, mix its versions in, and emit its
    per-contig classification output alongside mob_recon's.

    mob_recon's own script decompresses a gzipped assembly internally, unlike
    spaTyper/staphopia-sccmec in subworkflows/local/typing - no separate
    GUNZIP step is needed here.
*/

include { MOBSUITE_RECON } from '../../../modules/nf-core/mobsuite/recon/main'

workflow PLASMID_CLASSIFICATION {
    take:
    assembly // tuple val(meta), path(fasta), possibly gzipped

    main:
    ch_versions = channel.empty()

    MOBSUITE_RECON ( assembly )

    emit:
    mobsuite_contig_report    = MOBSUITE_RECON.out.contig_report    // channel: [ val(meta), path(contig_report.txt) ] - per-contig chromosome/plasmid-group call
    mobsuite_chromosome       = MOBSUITE_RECON.out.chromosome       // channel: [ val(meta), path(chromosome.fasta) ]
    mobsuite_plasmids         = MOBSUITE_RECON.out.plasmids         // channel: [ val(meta), path(plasmid_*.fasta) ]
    mobsuite_mobtyper_results = MOBSUITE_RECON.out.mobtyper_results // channel: [ val(meta), path(mobtyper_results.txt) ] - replicon/relaxase typing per plasmid, feeds stage 5b (#25)
    versions                  = ch_versions
}
