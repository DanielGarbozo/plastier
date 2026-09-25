# Stage 7 validation: methodology and results

Closed-genome benchmark of the four-tier framework (issues #22, #32, #33, #34, #36). The
simulated-read sensitivity analysis (#35) is separate and diagnostic only; it is described at the
end and reported in its own section once run. Raw result tables are in
[`stage7_results/`](stage7_results/); everything below regenerates with
`assets/validation/run_stage7.sbatch`.

**One-paragraph summary.** On 25 closed *S. aureus* genomes the tier framework called all 336 ARG
hits correctly (321 chromosomal, 14 High-confidence plasmid, 1 Moderate-confidence plasmid, no
Ambiguous). So did MOB-suite and RFPlasmid on their own; Platon alone missed one plasmid-borne
gene. **This benchmark cannot show that the framework beats a single tool**: it is at ceiling, the
positive class is only 15 hits from 8 genomes, and closed genomes make plasmid classification
easy. The evidence-integration logic is not contradicted, but it is also barely tested by this
data. PlasEval, which scores plasmid *binning*, is where the classifier imperfections show.

## Data

- **Ground truth:** `assets/validation/closed_genome_ground_truth.csv` (#31): 25 complete
  RefSeq *S. aureus* genomes and their 42 plasmid replicons, each with its accession and size.
  Every genome has at least one plasmid, so there are no chromosome-only controls.
- **Input to the pipeline:** the genome FASTA itself, via `--assemblies`
  (`assets/validation/fetch_reference_genomes.py`), so each ARG hit carries the replicon's own
  accession as its contig id. Stages 1-2 (reads, assembly) are skipped on purpose: a re-assembly
  renames contigs and the truth could no longer be looked up per contig.
- **Headers were rewritten** to Unicycler style (`>ACC length=N depth=1.00x circular=true`),
  sequences untouched. Stage 5c reads circularity and depth from exactly those fields and an
  NCBI header has neither, which would make the High-confidence tier impossible. These are
  complete genomes, so `circular=true` is a fact; it is set on **every** replicon, chromosome
  included, so it does not reveal which are plasmids. `depth=1.00x` gives no coverage evidence.
  **The High-confidence results depend on this choice** (`--keep-ncbi-headers` turns it off).

## Method

- **Pipeline:** stages 3-6, `-profile validation,singularity,slurm`. ARG calls from ABRicate,
  AMRFinderPlus and fARGene (RGI, DeepARG and argnorm are skipped, as in the `test` profile);
  MOB-suite, Platon and RFPlasmid; MLST, spa, SCCmec typing and SCCmecExtractor; then stage 5.
  753 tasks, none failed. AMRFinderPlus `--plus` is off, so metal/biocide genes that the ground
  truth lists for some plasmids (cadC, merA, qacA ...) are not among the calls.
- **Unit of evaluation:** one ARG *hit*, i.e. a distinct (gene, contig, start, stop) call from
  any ARG tool: 336 hits. A locus reported by several tools under different names counts once per
  name/coordinates, so heavily redundant genes weigh more.
- **Truth per hit** (`bin/evaluate_metrics.py`): on a plasmid if its contig accession is one of
  that genome's `Plasmid_Accession` values, otherwise chromosomal. This is per contig, not per
  gene name, so a gene carried on both a plasmid and the chromosome is scored correctly for each.
- **Metrics are per tier and never pooled.** For a plasmid tier the positive class is
  "truly plasmid", for Chromosomal it is "truly chromosome". FN counts true members of the class
  not assigned that tier, so recall is the share of *all* true plasmid hits that landed in this
  tier (which is why the Moderate row's recall is low: its 14 "misses" are in High).
  Ambiguous is an abstention with no true class, so only its call count is reported.
- **Baseline (#34):** each classifier's own `*_call` per hit, against the same truth. The tier
  framework is reduced to plasmid / chromosome / no call (Ambiguous). Coverage is reported so an
  abstention costs recall instead of hiding.
- **PlasEval (#33):** MOB-suite's plasmid bins (`primary_cluster_id`, qualified by genome)
  against the closed genomes' plasmids, with PlasEval `eval` (commit `8238356`).

## Results

### Per-tier precision / recall / F1 (#32)

| Tier | Calls | True plasmid | True chromosome | Precision | Recall | F1 |
| --- | --- | --- | --- | --- | --- | --- |
| High-confidence plasmid | 14 | 14 | 0 | 1.000 | 0.933 | 0.966 |
| Moderate-confidence plasmid | 1 | 1 | 0 | 1.000 | 0.067 | 0.125 |
| Chromosomal | 321 | 0 | 321 | 1.000 | 1.000 | 1.000 |
| Ambiguous | 0 | - | - | - | - | - |

No hit was placed in the wrong class. Of the 15 plasmid hits, 14 are High and 1 is Moderate
(Platon called that contig chromosomal, so agreement was 2 of 3). Nothing was Ambiguous, so the
abstention path was not exercised.

### Single-tool baseline vs the framework (#34)

| Method | Plasmid P / R | Chromosome P / R | Coverage |
| --- | --- | --- | --- |
| Tier framework | 1.000 / 1.000 | 1.000 / 1.000 | 1.00 |
| MOB-suite alone | 1.000 / 1.000 | 1.000 / 1.000 | 1.00 |
| RFPlasmid alone | 1.000 / 1.000 | 1.000 / 1.000 | 1.00 |
| Platon alone | 1.000 / 0.933 | 0.997 / 1.000 | 1.00 |

The framework matches the best single tools and is better than Platon alone by one hit. That is a
tie, not a win.

### PlasEval bin-level comparison (#33)

Overall (contig-level / base-pair-weighted): precision 0.875 / 0.899, recall 0.952 / 0.988,
F1 0.912 / 0.942. This is where MOB-suite is not perfect:

- three predicted bins each merge several true plasmids into one MOB-cluster (precision 0.33,
  0.33 and 0.50: `GCF_026625265.1|AA849`, `GCF_026625305.1|AB542`, `GCF_044789125.1|AC333`);
- two true plasmids (`NZ_CP170680.1`, `NZ_CP113006.1`) were not binned as plasmids at all
  (recall 0).

These errors are invisible in the per-tier table because none of those plasmids carried an ARG
that was misclassified. The tier framework does not repair bin-level errors; it only labels ARGs.

## What this does and does not show

- **Ceiling effect.** Every method scores near 1.0, so the benchmark has no power to rank them.
  Closed replicons are the easy case: each plasmid is one contig. The framework exists for
  fragmented draft assemblies, which is what the sensitivity analysis (#35) approximates.
- **Small positive class.** 15 plasmid hits, 8 genomes, 8 distinct genes (erm(C) alone is 5).
  One hit changes plasmid recall by about 7 points. No confidence intervals are given because
  the counts are too small to make them meaningful.
- **Circularity was supplied**, not detected (see Data). Without it no hit can be High.
- **Ambiguous and the length floor were not exercised.** No hit was Ambiguous, and no ARG-bearing contig
  is under the 1,000 bp floor (the smallest is 2,366 bp). Those rules are untested here.
- **No rule change is proposed.** The coverage-ratio (1.5x) and contig-length-floor (1,000 bp)
  defaults in `docs/decisions.md` are unchanged; this data neither supports nor argues against
  them, since the coverage ratio was fixed at 1.0x and no ARG-bearing contig is short.
- **Single species, single database version.** Nothing here speaks to other species.

## Sensitivity analysis (#35): diagnostic only

`assets/validation/run_sensitivity.sbatch` simulates paired 2x150 reads (ART, HiSeq 2500 profile)
from the genomes with the most plasmid-borne ARGs, at several depths, base-quality shifts and
plasmid copy numbers, assembles them, and compares each gene's tier with its closed-genome tier
(`bin/sensitivity_summary.py`). Contigs are renamed by assembly, so genes are matched by name.
Plasmid copy number is an assumption of the simulation (p1 = same depth as the chromosome, so
stage 5's coverage evidence is absent; p3 = a plausible multi-copy plasmid). **This is never an
accuracy figure**; it shows how much depth and quality move the calls.

## Reproducing

```bash
sbatch --account=<acct> --partition=<partition> assets/validation/run_stage7.sbatch   # run + score
assets/validation/score_stage7.sh <OUTDIR>                                            # re-score only
pytest tests/python                                                                   # the scoring code
```
