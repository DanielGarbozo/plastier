# plastier roadmap

Tracks phase scope and approval status for plastier as a project —
distinct from `CHANGELOG.md`, which tracks released pipeline versions.

## Phase 0 — Pilot (current)

- **Scope:** 10 *S. aureus* genomes, stages 1–2 only (nf-core/fetchngs →
  nf-core/bacass).
- **Status:** Approved.
- **Execution:** AWS execution on the paid codespace is restricted to
  Nkiruka; see `CLAUDE.md`.
- **Local/CI scope:** `-profile test,docker` only — near-zero cost, runs
  anywhere.
- **HPC scope:** `-profile test,singularity,slurm` via
  `assets/slurm/smoke_test.sbatch`. The smoke test passes end to end
  (stages 1–6 on the `test` profile, 2026-09-25).

## Phase 1 — Plasmid classification + evidence integration (implemented)

- **Scope:** MOB-suite, Platon, RFPlasmid (stage 4) and the local
  evidence-integration subworkflow that resolves ARG-to-contig calls into
  the four-tier framework (stage 5).
- **Status:** Approved 2026-09-25. Stages 4–5 are implemented and wired into
  `workflows/plastier.nf`.

## Phase 2 — Typing + validation (in progress)

- **Scope:** MLST/*spa*/SCCmec typing (stage 6), closed-genome benchmark,
  PlasEval comparison, single-tool baseline comparator, simulated-read
  sensitivity analysis (stage 7).
- **Status:** Approved 2026-09-25.
  - Stage 6 typing: implemented.
  - Stage 7 run 2026-09-27 on 25 closed genomes (#31-#36, PR #72): per-tier
    precision/recall/F1 (#32), single-tool baseline (#34) and PlasEval
    bin-level comparison (#33) all done and documented in
    `docs/stage7_validation.md`. Every method scores at or near 1.0 - the
    benchmark is at ceiling (only 15 plasmid hits, 8 genomes) and cannot show
    the tier framework beats a single tool; the doc says so explicitly. The
    simulated-read sensitivity analysis (#35, diagnostic only) found the
    plasmid tier is unstable when the simulated plasmid has no coverage
    advantage over the chromosome, one false plasmid call at 10x depth, and
    outright assembly failure at the harshest quality setting tested.

## Phase 3 — Scale-up (in progress)

- **Scope:** Move beyond the 10-genome pilot to a discovery run of about 20
  *S. aureus* genomes, executed on the Slurm cluster (`-profile slurm`).
- **Status:** Approved in principle 2026-09-25 for ~20 genomes. Candidates
  selected (`docs/phase3_candidates.md`) and a full run completed 2026-09-27
  (PR #72); per-sample ST/*spa*/SCCmec/ARG-tier results are in
  `docs/phase3_results.md` for review before `conf/full.config`'s `sra_ids`
  is set (still unset - fails closed). One sample has no MLST ST and needs a
  closer look; SCCmecExtractor only resolved the physical SCCmec cassette in
  2/20 samples (fragmented drafts split the *att* sites across contigs), so
  the stage 5e SCCmec-override path got little exercise on this cohort.

## Out of scope until explicitly requested

- Any cloud execution beyond what Nkiruka runs on the paid codespace.
- Any species other than *S. aureus* for the pilot (the workflow is
  designed to be species-agnostic, but this has not been validated yet).

## Change log for this document

- 2026-07-31: Initial stub created alongside stages 1–2 scaffold.
- 2026-09-25: Dropped the Phase 0 budget cap; brought phase statuses in line
  with the code (stages 4–6 implemented); Phases 1–3 approved to proceed,
  Phase 3 at ~20 genomes on Slurm.
- 2026-09-27: Stage 7 validation run and Phase 3 scale-up run both completed
  (PR #72); results and open review questions recorded in
  `docs/stage7_validation.md` and `docs/phase3_results.md`.
