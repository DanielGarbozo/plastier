/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC                   } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { paramsSummaryMultiqc      } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML    } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText    } from '../subworkflows/local/utils_nfcore_plastier_pipeline'
include { validateInputSamplesheet  } from '../subworkflows/local/utils_nfcore_plastier_pipeline'
include { FETCHNGS                  } from '../subworkflows/local/fetchngs/main'
include { BACASS                    } from '../subworkflows/local/bacass/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PLASTIER {

    take:
    ch_samplesheet // channel: samplesheet read in from --input (may be empty if --sra_ids is set)
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    //
    // STAGE 1 (optional): retrieve public reads via nf-core/fetchngs (--sra_ids)
    //
    // If --sra_ids is given, download the runs it lists and build the stage-2
    // input samplesheet from the result, in place of a pre-existing --input file.
    // Fails fast if neither is provided, since bacass has nothing to assemble.
    //
    def ch_bacass_input
    if (params.sra_ids) {
        def ch_ids = channel
            .fromPath(params.sra_ids, checkIfExists: true)
            .splitCsv(header: false, sep: '', strip: true)
            .map { it[0] }
            .unique()

        FETCHNGS ( ch_ids )
        ch_versions = ch_versions.mix(FETCHNGS.out.versions)

        ch_bacass_input = FETCHNGS.out.samplesheet
            .flatMap { csv -> samplesheetToList(csv, "${projectDir}/assets/schema_input.json") }
            .map { meta, fastq_1, fastq_2 ->
                fastq_2
                    ? [ meta.id, meta + [ single_end: false ], [ fastq_1, fastq_2 ] ]
                    : [ meta.id, meta + [ single_end: true ], [ fastq_1 ] ]
            }
            .groupTuple()
            .map { samplesheet -> validateInputSamplesheet(samplesheet) }
            .map { meta, fastqs -> [ meta, fastqs.flatten() ] }
    } else {
        if (!params.input) {
            error("Either --input (a fastq samplesheet) or --sra_ids (a list of public accessions for nf-core/fetchngs) must be provided.")
        }
        ch_bacass_input = ch_samplesheet
    }

    //
    // STAGE 2: uniform QC and assembly via nf-core/bacass
    //
    BACASS ( ch_bacass_input )
    ch_versions = ch_versions.mix(BACASS.out.versions)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name:  'plastier_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'plastier'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )
    // MULTIQC.out.report: pipeline-level summary (params, versions, methods description)
    // BACASS.out.multiqc_report: stage-2 per-sample QC/assembly report
    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> report }
        .mix(BACASS.out.multiqc_report.flatten())
        .toList() // channel: [ /path/to/multiqc_report.html, ... ]
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
