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

**plastier/plastier** is a bioinformatics pipeline that assigns every detected resistance gene to a plasmid or to the chromosome, states a confidence level for that assignment, and reports how accurate each confidence level actually is when tested against genomes where the answer is known. 

It processes raw sequencing data through five core stages:

1. **Data Retrieval (Stage 1):** Downloads raw short-reads from public databases via `nf-core/fetchngs` (optional).
2. **Assembly & Annotation (Stage 2):** Performs *de novo* genome assembly and structural annotation using `nf-core/bacass` (Unicycler/SPAdes, Prokka).
3. **AMR Detection (Stage 3):** Scans the assemblies for resistance genes using an arsenal of tools via `funcscan` (AMRFinderPlus, FARGENE, RGI, DeepARG, Abricate) and standardizes them using hAMRonization.
4. **Plasmid Classification & Typing (Stage 4):** Evaluates contig origin using multiple classifiers (Platon, MOB-suite, RFPlasmid) and performs strain typing (MLST, spaTyper, SCCmec).
5. **Evidence Integration (Stage 5):** A custom rules engine (`tier_resolution.py`) evaluates the conflicting outputs of the classifiers and assigns a final 4-tier confidence level to every detected ARG (High-confidence plasmid, Moderate-confidence plasmid, Chromosomal, or Ambiguous).


## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/get_started/environment_setup/overview) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/get_started/run-your-first-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,fastq_1,fastq_2
isolate_1,reads/isolate_1_R1.fastq.gz,reads/isolate_1_R2.fastq.gz
isolate_2,reads/isolate_2_R1.fastq.gz,reads/isolate_2_R2.fastq.gz
```

Each row represents a pair of fastq files (paired-end). Alternatively, if you are starting from public data, you can provide an SRA IDs file (CSV with an `id` column) to automatically download the reads using the `fetchngs` module.

Now, you can run the pipeline using:

```bash
nextflow run plastier/plastier \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/running/run-pipelines#using-parameter-files).

## Credits

plastier/plastier was originally written by Daniel Garbozo.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](docs/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use plastier/plastier for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
