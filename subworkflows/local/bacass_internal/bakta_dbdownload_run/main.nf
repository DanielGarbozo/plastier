//
// Annotation of Bacterial genomes with Bakta
//

include { BAKTA_BAKTADBDOWNLOAD } from '../../../../modules/nf-core/bakta/baktadbdownload'
include { UNTAR                 } from '../../../../modules/nf-core/untar'
include { BAKTA_BAKTA           } from '../../../../modules/nf-core/bakta/bakta'


workflow BAKTA_DBDOWNLOAD_RUN {
    take:
    ch_fasta                // channel: [ val(meta), path(fasta)  ]
    ch_path_baktadb         // channel: [ path(fasta) ]
    val_baktadb_download    // value: boolean

    main:
    ch_versions = channel.empty()

    //
    // SUBWORKFLOW: Parse, download and/or untar Bakta database
    //
    if( ch_path_baktadb ){
        if (ch_path_baktadb.endsWith('.tar.gz')){
            ch_baktadb_tar  = channel.from(ch_path_baktadb).map{ db -> [ [id: 'baktadb'], db ]}

            // MODULE: untar database
            // NOTE: deviation from upstream bacass 2.6.1 - the vendored modules/nf-core/untar
            // and modules/nf-core/bakta/* emit software versions via the `versions` topic
            // instead of a classic `.out.versions` channel, already auto-collected by
            // BACASS's topic-versions block, so no explicit mix() is needed here anymore.
            UNTAR( ch_baktadb_tar )
            ch_path_baktadb = UNTAR.out.untar.map{ _meta, db -> db }
        }
    } else if (!ch_path_baktadb && val_baktadb_download){
        // MODULE: Downlado Bakta database from zenodo
        BAKTA_BAKTADBDOWNLOAD()
        ch_path_baktadb  = BAKTA_BAKTADBDOWNLOAD.out.db

    } else if (!ch_path_baktadb && !val_baktadb_download ){
        exit 1, "The Bakta database argument is missing. To enable the workflow to access the Bakta database, please include the path using '--baktadb' or use '--bakdtadb_download true' to download the Bakta database."
    }

    //
    // MODULE: BAKTA, gene annotation
    //
    // Setup input channel for Bakta process
    BAKTA_BAKTA (
        ch_fasta,
        ch_path_baktadb,
        [],
        [],
        [],
        []
    )
    ch_bakta_txt_multiqc    = BAKTA_BAKTA.out.txt

    emit:
    versions                = ch_versions.ifEmpty(null) // channel: [ path(versions.yml) ]
    bakta_txt_multiqc       = ch_bakta_txt_multiqc      // channel: [ meta, path(*.txt)  ]
}
