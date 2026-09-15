# plastier/plastier: Output

## Introduction

This document describes the output produced by the pipeline. Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

The directories listed below will be created in the results directory (`--outdir`) after the pipeline has finished. All paths are relative to the top-level results directory. Each tool publishes to a directory named after itself (e.g. `unicycler/`, `mobsuite/`), so the layout below is grouped by pipeline stage rather than by literal directory nesting.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:

- [Stage 1: Read retrieval](#stage-1-read-retrieval-optional) (optional) - download public runs via nf-core/fetchngs
- [Stage 2: QC and assembly](#stage-2-qc-and-assembly) - read QC, trimming, assembly, and annotation via nf-core/bacass
- [Stage 3: ARG screening](#stage-3-arg-amr-gene-screening) - antimicrobial resistance gene detection via nf-core/funcscan
- [Stage 4: Plasmid classification](#stage-4-plasmid-classification) - per-contig plasmid/chromosome calls from MOB-suite, Platon, and RFPlasmid
- [Stage 5: Evidence integration](#stage-5-evidence-integration) - the four-tier confidence call per ARG, and the intermediate signals behind it
- [Stage 6: Strain typing](#stage-6-strain-typing) - MLST, *spa* typing, SCCmec typing, and cassette extraction
- [MultiQC](#multiqc) - Aggregate report describing results and QC from the whole pipeline
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

> [!NOTE]
> Stage 7 (validation/benchmarking) is not yet implemented (see [issue #22](https://github.com/DanielGarbozo/plastier/issues/22)) and produces no output yet.

## Stage 1: Read retrieval (optional)

Only runs when `--sra_ids` is supplied instead of `--input`. Wraps [nf-core/fetchngs](https://nf-co.re/fetchngs) to download public runs and build the samplesheet that stage 2 consumes.

<details markdown="1">
<summary>Output files</summary>

- `fetchngs/` (or the equivalent tool-named subdirectories, e.g. `sra/`)
  - `samplesheet/samplesheet.csv`: samplesheet built from the downloaded runs, in nf-core/fetchngs's own `sample,fastq_1,fastq_2,...` format - used internally to build stage 2's input, not the same schema as `assets/schema_input.json`.
  - `*.fastq.gz`: downloaded raw reads.
  - `id_mappings.csv` / `multiqc_config.yml`: ENA/SRA metadata mappings, when `--sample_mapping_fields` is set.

</details>

## Stage 2: QC and assembly

Wraps [nf-core/bacass](https://nf-co.re/bacass) for read QC/trimming, genome assembly, and annotation. Currently exercised in short-read mode (`--assembly_type short`).

<details markdown="1">
<summary>Output files</summary>

- `fastqc/`
  - `*_fastqc.html` / `*_fastqc.zip`: raw-read quality metrics (see [FastQC](#fastqc) below).
- `fastp/`
  - `*.fastp.json` / `*.fastp.html`: adapter/quality trimming report per sample.
- `unicycler/` (or the configured assembler's own directory, e.g. `canu/`, `flye/`)
  - `*.assembly.fasta.gz`: the assembled genome (input to every downstream stage).
  - `*.assembly.gfa.gz`: assembly graph.
- `quast/`
  - Assembly quality metrics (N50, contig counts, genome length).
- `prokka/` (or `bakta/`/`dfast/`, depending on `--annotation_tool`)
  - `*.faa`: predicted protein sequences - the direct input to stage 3 (ARG screening).
  - `*.gff`: genome annotation.
- `kraken2/`, `kmerfinder/` (unless skipped): contamination/species-identity screening.
- `busco/` (unless skipped): genome completeness assessment.

</details>

The assembly (`*.assembly.fasta.gz`) and annotation (`*.faa`/`*.gff`) produced here are what every subsequent stage operates on - stages 3-6 never re-annotate or re-assemble.

## Stage 3: ARG (AMR gene) screening

Wraps [nf-core/funcscan](https://nf-co.re/funcscan)'s ARG module, running AMRFinderPlus, RGI, ABRicate, DeepARG, and fARGene against the stage-2 annotation (each individually skippable via `--arg_skip_<tool>`).

<details markdown="1">
<summary>Output files</summary>

- `amrfinderplus/`, `rgi/`, `abricate/`, `deeparg/`, `fargene/`: per-tool raw resistance-gene calls.
- `hamronization/`
  - `hamronization_combined_report.tsv`: all tools' calls normalised into one schema by [hAMRonization](https://github.com/pha4ge/hAMRonization) - this is the report stage 5 (evidence integration) reads `gene_symbol`/`input_sequence_id` from to key every downstream signal.
- `argnorm/` (unless `--arg_skip_argnorm`): resistance gene nomenclature normalised against the CARD ontology.

</details>

## Stage 4: Plasmid classification

Runs MOB-suite, Platon, and RFPlasmid on the stage-2 assembly to independently call each contig as plasmid or chromosome. All three are individually skippable (`--plasmid_skip_platon`, `--plasmid_skip_rfplasmid`; MOB-suite always runs).

<details markdown="1">
<summary>Output files</summary>

- `mobsuite/`
  - `*_contig_report.txt`: per-contig plasmid/chromosome call and plasmid cluster grouping.
  - `*_mobtyper_results.txt`: replicon and relaxase typing per predicted plasmid - feeds the mobility-marker signal in stage 5.
  - `*_chromosome.fasta` / `*_plasmid_*.fasta`: contigs split out by MOB-suite's own call.
- `platon/` (unless skipped)
  - `*.tsv`: per-contig plasmid/chromosome call with supporting evidence.
  - `*.json`: full per-contig detail.
- `rfplasmid/` (unless skipped)
  - `prediction.csv`: per-contig plasmid/chromosome call (random-forest classifier).
  - `prediction_full.csv`: full classifier output including per-class probabilities.

</details>

## Stage 5: Evidence integration

Combines the stage 3 ARG report with the stage 4 classifier calls, stage 6 typing, and two structural checks (circularity/coverage, contig length) into the pipeline's core deliverable: one confidence tier per detected resistance gene. See [`docs/decisions.md`](decisions.md) for the full rule table and the reasoning behind each threshold.

<details markdown="1">
<summary>Output files</summary>

- `classifieragreement/`
  - `*.classifier_agreement.tsv`: whether MOB-suite, Platon, and RFPlasmid agree on each ARG's contig (`unanimous_plasmid`, `majority_plasmid`, `unanimous_chromosome`, `majority_chromosome`, or `disagreement`), plus each classifier's individual call and vote counts.
- `mobilitymarkers/`
  - `*.mobility_markers.tsv`: replicon/relaxase type and predicted mobility from MOB-suite's typer output, and whether any marker (`has_marker`) was found for the contig.
- `circularitycoverage/`
  - `*.circularity_coverage.tsv`: whether the contig is circularised, its coverage ratio relative to the assembly median, and whether that ratio matches plasmid-like coverage (`--plasmid_coverage_ratio_threshold`, default 1.5x).
- `contiglengthfloor/`
  - `*.contig_length_floor.tsv`: contig length and whether it falls below the trusted-call floor (`--plasmid_min_contig_length`).
- `sccmecoverride/`
  - `*.sccmec_override.tsv`: whether the ARG sits inside a typed SCCmec cassette, and the cassette span if so.
- `tierresolution/`
  - **`*.tier_resolution.tsv`: the final output.** One row per detected ARG, keyed by `gene_symbol` and `input_sequence_id`, with every upstream signal above plus the resolved `tier` column - one of `High-confidence plasmid`, `Moderate-confidence plasmid`, `Ambiguous`, or `Chromosomal`.

</details>

## Stage 6: Strain typing

Runs MLST, *spa* typing (SPAtyper), SCCmec typing (StaphopiaSCCmec), and SCCmec cassette extraction on the stage-2 assembly. Separates horizontal AMR transfer (stages 3-5) from clonal strain background.

<details markdown="1">
<summary>Output files</summary>

- `mlst/`
  - `*.tsv`: multi-locus sequence type call.
- `spatyper/`
  - `*.tsv`: *spa* type call (*S. aureus*-specific).
- `staphopiasccmec/`
  - `*.tsv`: SCCmec type call (*S. aureus*-specific).
- `sccmecextractor/`
  - `sccmec_unified_report.tsv`: coordinates of any typed SCCmec cassette on the assembly - consumed directly by stage 5's SCCmec override.

</details>

### MultiQC

<details markdown="1">
<summary>Output files</summary>

- `multiqc/`
  - `multiqc_report.html`: a standalone HTML file that can be viewed in your web browser.
  - `multiqc_data/`: directory containing parsed statistics from the different tools used in the pipeline.
  - `multiqc_plots/`: directory containing static images from the report in various formats.

</details>

[MultiQC](http://multiqc.info) is a visualization tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in the report data directory.

Results generated by MultiQC collate pipeline QC from supported tools e.g. FastQC. The pipeline has special steps which also allow the software versions to be reported in the MultiQC output for future traceability. For more information about how to use MultiQC reports, see <http://multiqc.info>.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.

</details>

[Nextflow](https://docs.seqera.io/platform-cloud/reports/overview) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
