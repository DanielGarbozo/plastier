# plastier roadmap

Tracks phase scope, approval status, and budget for plastier as a project —
distinct from `CHANGELOG.md`, which tracks released pipeline versions.

## Phase 0 — Pilot (current)

- **Scope:** 10 *S. aureus* genomes, stages 1–2 only (nf-core/fetchngs →
  nf-core/bacass).
- **Status:** Approved.
- **Budget:** $8 USD AWS spend, capped. Execution on the paid AWS codespace
  is restricted to Nkiruka; see `CLAUDE.md`.
- **Local/CI scope:** `-profile test,docker` only — near-zero cost, runs
  anywhere.

## Phase 1 — Plasmid classification + evidence integration (not started)

- **Scope:** Add MOB-suite, Platon, RFPlasmid (stage 4) and the local
  evidence-integration subworkflow that resolves ARG-to-contig calls into
  the four-tier framework (stage 5).
- **Status:** Not approved. Requires sign-off before implementation begins.

## Phase 2 — Typing + validation (not started)

- **Scope:** MLST/*spa*/SCCmec typing (stage 6), closed-genome benchmark,
  PlasEval comparison, single-tool baseline comparator, simulated-read
  sensitivity analysis (stage 7).
- **Status:** Not approved.

## Phase 3 — Scale-up (not started)

- **Scope:** Move beyond the 10-genome pilot to a full discovery run.
- **Status:** Not approved. Requires an explicit budget decision separate
  from the Phase 0 cap — do not assume Phase 0's $8 cap extends to this
  phase.

## Out of scope until explicitly requested

- Any cloud execution beyond what Nkiruka runs on the paid codespace.
- Any species other than *S. aureus* for the pilot (the workflow is
  designed to be species-agnostic, but this has not been validated yet).

## Change log for this document

- 2026-07-31: Initial stub created alongside stages 1–2 scaffold.
