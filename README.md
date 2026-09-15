# plastier/plastier

[![GitHub Actions CI Status](https://github.com/plastier/plastier/actions/workflows/nf-test.yml/badge.svg)](https://github.com/plastier/plastier/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/plastier/plastier/actions/workflows/linting.yml/badge.svg)](https://github.com/plastier/plastier/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A525.10.4-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-4.1.0-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/4.1.0)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/plastier/plastier)

## Introduction

**plastier/plastier** is a bioinformatics pipeline that assigns every detected antimicrobial resistance (AMR) gene to a plasmid or to the chromosome, states a confidence level for that assignment (a four-tier framework), and reports how accurate each confidence level actually is when tested against genomes where the answer is known. Distinguishing plasmid-borne from chromosomal resistance genes matters because plasmid-borne AMR spreads horizontally between strains and species, while chromosomal AMR is inherited clonally - the two have very different implications for surveillance and outbreak response.

The pipeline takes raw short-read (Illumina) sequencing data - either supplied directly or fetched from public archives - through assembly, AMR gene screening, plasmid classification, strain typing, and evidence integration:

1. **Read retrieval** (optional, `--sra_ids`) - download public runs from SRA/ENA/DDBJ/GEO ([nf-core/fetchngs](https://nf-co.re/fetchngs))
2. **QC and assembly** - read QC, trimming, and genome assembly ([nf-core/bacass](https://nf-co.re/bacass))
3. **AMR gene (ARG) screening** - detect resistance genes with AMRFinderPlus, RGI, ABRicate, DeepARG, and fARGene ([nf-core/funcscan](https://nf-co.re/funcscan))
4. **Plasmid classification** - classify assembly contigs as plasmid or chromosome with MOB-suite, Platon, and RFPlasmid
5. **Strain typing** - MLST, *spa* typing, SCCmec typing, and SCCmec cassette extraction
6. **Evidence integration** - resolve each ARG-to-contig call into the four-tier confidence framework by combining signals from steps 3-5
7. **Validation** (planned) - benchmark tier accuracy against closed genomes with known plasmid content

Stages 1-6 are implemented; stage 7 (validation) is tracked in [issue #22](https://github.com/DanielGarbozo/plastier/issues/22) and its sub-issues. See [`docs/roadmap.md`](docs/roadmap.md) for current project phase and approval status, and [`docs/decisions.md`](docs/decisions.md) for the reasoning behind each tier-resolution rule.

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data. `plastier` consumes [nf-core/bacass](https://nf-co.re/bacass)'s own samplesheet schema, since stage 2 is the immediate consumer of `--input`:

`samplesheet.csv`:

```csv
ID,R1,R2,LongFastQ,Fast5,GenomeSize
SAMPLE_1,/path/to/fastq/SAMPLE_1_R1.fastq.gz,/path/to/fastq/SAMPLE_1_R2.fastq.gz,NA,NA,2.8m
SAMPLE_2,/path/to/fastq/SAMPLE_2_R1.fastq.gz,/path/to/fastq/SAMPLE_2_R2.fastq.gz,NA,NA,2.8m
```

Each row represents one sample's paired-end short reads. `LongFastQ`, `Fast5`, and `GenomeSize` may be left as `NA` for short-read-only assembly (the pipeline's currently supported mode); see [`assets/schema_input.json`](assets/schema_input.json) for the full column specification and [`docs/usage.md`](docs/usage.md) for details.

Alternatively, in place of `--input`, you can supply `--sra_ids <path/to/ids.txt>` - a plain text file of public SRA/ENA/DDBJ/GEO accessions, one per line - and the pipeline will fetch and assemble those runs directly, skipping the need to build a samplesheet by hand.

Now, you can run the pipeline using:

```bash
nextflow run plastier/plastier \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

For more details and full parameter documentation, see [`docs/usage.md`](docs/usage.md) and [`docs/output.md`](docs/output.md) for a description of the pipeline's output files.

## Credits

plastier/plastier was originally written by [Daniel Garbozo](https://github.com/DanielGarbozo).

This pipeline builds on top of, and integrates, several existing nf-core pipelines and community-maintained tools; see the [Citations](#citations) section below for the full list.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md). Open issues are tracked on the [Issues board](https://github.com/DanielGarbozo/plastier/issues); [`docs/roadmap.md`](docs/roadmap.md) describes which project phases are currently approved for implementation.

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
