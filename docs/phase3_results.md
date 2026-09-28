# Phase 3 scale-up: results

The 20 candidate runs from [`phase3_candidates.md`](phase3_candidates.md) run end-to-end
(`assets/phase3/run_phase3.sbatch`, Slurm, `-profile full,singularity,slurm`). Per-sample
outputs are under `<outdir>/{mlst,spatyper,staphopiasccmec,sccmecextractor,tier}`; the summary
table below is [`phase3_results/phase3_summary.tsv`](phase3_results/phase3_summary.tsv), built
by pulling MLST species/ST, *spa* type, staphopia-sccmec's SCCmec call, SCCmecExtractor's
cassette-resolution status, and the tier-resolution ARG counts for each sample.

**This is a review checkpoint, not a validated result set** - see "Open questions" below before
treating any of it as final or wiring `sra_ids` into `conf/full.config`.

## Summary

| Sample | ST | spa | SCCmec (staphopia) | Cassette resolved? | ARG (High/Mod/Amb/Chrom) |
|---|---|---|---|---|---|
| DRR456378 | 15 | t084 | - | no (no_ccr) | 5 (3/0/0/2) |
| DRR484144 | 380 | t12889 | IV+IVc | no (cross_contig) | 7 (0/3/1/3) |
| DRR484283 | 5 | t002 | II+IIa | no (cross_contig) | 12 (3/0/6/3) |
| ERR7434454 | 398 | t011 | V+VII | no (cross_contig) | 7 (0/1/1/5) |
| SRR12333138 | 97 | t267 | mecA only | no (cross_contig) | 8 (4/0/2/2) |
| SRR13450679 | 1 | t398 | - | no (no_ccr) | 1 (0/0/0/1) |
| SRR14267647 | 5 | t002 | mecA only | yes | 10 (4/0/0/6) |
| SRR14862042 | 5 | t002 | IV+IVc | no (cross_contig) | 10 (0/0/7/3) |
| SRR18339536 | 944 | 04-34-17-23-24-24-17 | IIa | yes | 4 (0/0/0/4) |
| SRR2057031 | 93 | t202 | IV+IVa | no (cross_contig) | 5 (3/0/0/2) |
| SRR21098763 | 1 | t127 | - | no (no_ccr) | 4 (3/0/0/1) |
| SRR23309788 | 59 | t216 | - | no (no_ccr) | 5 (1/0/0/4) |
| SRR26374374 | 1 | t127 | - | no (no_ccr) | 5 (4/0/0/1) |
| SRR27322154 | 30 | t012 | - | no (no_ccr) | 6 (0/0/1/5) |
| SRR29758699 | 97 | t267 | - | no (no_ccr) | 4 (0/0/0/4) |
| SRR36881360 | 5 | t002 | II+IIa | no (cross_contig) | 13 (3/0/2/8) |
| SRR3731507 | 1207 | t1778 | - | no (no_ccr) | 6 (3/0/0/3) |
| SRR6332104 | 1 | t127 | - | no (no_ccr) | 4 (3/0/0/1) |
| SRR6347285 | - | t878 | - | no (no_ccr) | 4 (3/0/0/1) |
| SRR7867492 | 239 | t030 | III+V+VII+IIIa | no (cross_contig) | 11 (0/0/2/9) |

ARG counts are ARG x contig-location rows from `tier_resolution.tsv` (a gene present on several
contigs/copies counts once per row), split by tier: **High**-confidence plasmid, **Mod**erate,
**Amb**iguous, **Chrom**osomal.

## Open questions for review

- **SRR6347285 has no ST** (all other 19 got one). MLST still calls the species *S. aureus* and
  it typed a *spa* (t878), so this reads as a novel/incomplete allele profile rather than
  necessarily a bad accession - worth a manual look (e.g. contamination, low coverage) before
  deciding whether to keep or swap it (`--exclude` in `select_candidates.py`, seed 42).
- **SCCmecExtractor resolved the physical cassette in only 2/20 samples** (SRR14267647,
  SRR18339536); the other 18 failed, split `no_ccr` (10 - no *ccr* recombinase genes found,
  either genuinely non-MRSA or the region wasn't assembled) and `cross_contig` (8 - the *attR*/
  *attL* sites landed on different contigs of a fragmented draft, so the element's coordinates
  can't be resolved). This means the stage 5e SCCmec-override branch of the tier framework
  (docs/decisions.md #28) essentially did not fire on this cohort - it needs assemblies where
  the SCCmec cassette isn't split across contigs, which short-read WGS drafts often don't give.
  This is expected of short-read assemblies and is not itself a bug, but it means Phase 3 as run
  gives little signal on that specific code path.
- **No accuracy claim.** There is no ground truth for these public isolates (unlike the closed
  Stage 7 genomes), so this table is for manual plausibility review (do the ST/spa/SCCmec/ARG
  combinations look like real *S. aureus* clones, e.g. ST5-t002 is a well-known
  hospital-associated lineage), not a benchmark result.

## Reproducing

```bash
sbatch --account=<acct> --partition=<partition> assets/phase3/run_phase3.sbatch
```
