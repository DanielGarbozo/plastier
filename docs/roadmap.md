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
  - Stage 7: ground-truth set (#31) done; per-tier metrics script (#32) and
    PlasEval converter (#33) exist but are not yet wired into the pipeline
    or run against real output. #34, #35 and #36 not started.

## Phase 3 — Scale-up (not started)

- **Scope:** Move beyond the 10-genome pilot to a discovery run of about 20
  *S. aureus* genomes, executed on the Slurm cluster (`-profile slurm`).
- **Status:** Approved in principle 2026-09-25 for ~20 genomes. The accession
  list still has to be curated and reviewed, and `conf/full.config` stays a
  stub until then.

## Out of scope until explicitly requested

- Any cloud execution beyond what Nkiruka runs on the paid codespace.
- Any species other than *S. aureus* for the pilot (the workflow is
  designed to be species-agnostic, but this has not been validated yet).

## Change log for this document

- 2026-07-31: Initial stub created alongside stages 1–2 scaffold.
- 2026-09-25: Dropped the Phase 0 budget cap; brought phase statuses in line
  with the code (stages 4–6 implemented); Phases 1–3 approved to proceed,
  Phase 3 at ~20 genomes on Slurm.
