# Decision log

Records changes to the four-tier evidence framework (see `CLAUDE.md`) and
other decisions that affect what the benchmark in stage 7 validates. Every
change to the tier rule logic needs an entry here explaining *why*, not just
a code diff - silently changing the rules invalidates prior validation
results.

Entries are appended chronologically, newest last. Do not edit or delete
past entries; if a decision is reversed, add a new entry that supersedes it
and say so explicitly.

## Template for new entries

```
## YYYY-MM-DD: Short title

**Decision:**
**Why:**
**Affects:** (which tier(s), which validation results, if any)
**Author:**
```

## 2026-07-31: Project scaffold created

**Decision:** Initialised the plastier repository as an nf-core-style
pipeline. Built stages 1–2 only (nf-core/fetchngs → nf-core/bacass), pinned
to fetchngs v1.12.0 and bacass v2.6.1. No changes to the four-tier framework
itself - stages 4-5 (plasmid classification, evidence integration) are not
yet implemented, so this entry exists only to mark the starting point.

**Why:** Establish a working, tested foundation (`-profile test,docker`)
before adding the plasmid-classification and tier-resolution logic that the
framework actually depends on.

**Affects:** None yet - no tier logic exists in the codebase at this point.

**Author:** Daniel Garbozo

## 2026-08-25: Coverage-ratio threshold for stage 5c (issue #26)

**Decision:** The "High-confidence plasmid" tier's coverage-ratio condition
("coverage ratio matches plasmid copy number") is implemented as: the
contig's Unicycler-reported depth (relative to the assembly's own median,
which Unicycler already normalises to ~1.00x for the chromosome) is >= 1.5x.
Configurable via `--plasmid_coverage_ratio_threshold` (script default 1.5,
wired into the pipeline as a real param - see
subworkflows/local/evidence_integration), not hardcoded.

**Why:** Real chromosomal contigs from this repo's own bacass test fixtures
vary roughly 0.97x-1.08x (see bin/circularity_coverage.py's docstring for the
actual header values this was checked against). 1.5x is a conservative floor
above that normal noise band - not a literal claim about true plasmid copy
number, which this pipeline has no direct way to measure. Chosen without
validation data (stage 7, issue #22, doesn't exist yet), so this is a
placeholder default, not a calibrated one.

**Affects:** High-confidence plasmid tier only (the circularisation branch of
the same OR condition is unaffected). Must be revisited once stage 7 has
real closed-genome ground truth to check this threshold against - if it
turns out too strict/lenient, that's a future entry here, not a silent
tweak.

**Author:** Claude (session with Daniel Garbozo)

## 2026-08-25: Contig-length floor default for stage 5d (issue #27)

**Decision:** Contigs shorter than 1000bp (Unicycler's `length=` header field)
are routed to the Ambiguous tier regardless of what the classifiers say.
Configurable via `--plasmid_min_contig_length`, not hardcoded.

**Why:** 1000bp is a commonly used floor in short-read plasmid-classification
literature below which assembly fragmentation and misassembly make a
plasmid/chromosome call unreliable. Not calibrated against this pipeline's
own data - no validation set exists yet (stage 7, issue #22).

**Affects:** Ambiguous tier only. Revisit once stage 7 has real data.

**Author:** Claude (session with Daniel Garbozo)

## 2026-08-25: SCCmecExtractor added for stage 5e (issue #28), not staphopia-sccmec

**Decision:** The SCCmec-cassette override needs the physical coordinates of
the cassette on its contig, to test whether an ARG's own coordinates fall
inside it. staphopia-sccmec (already in `subworkflows/local/typing` for
stage 6 typing) was checked against its real `--json` output and confirmed
to report only a genome-level presence/absence call per SCCmec type - no
contig, no position, nothing to overlap against. A second tool,
SCCmecExtractor (pip: `sccmecextractor`, pinned Docker image
`alisonmacfadyen/sccmecextractor:v1.6.0`), was added specifically for this -
it locates the actual att-site attachment boundaries per contig. Verified by
running it directly (not just reading its docs) against its own bundled test
genome before committing to it.

Both tools stay in stage 6: staphopia-sccmec for typing (already relied on
elsewhere), SCCmecExtractor for the coordinates stage 5e needs. This is a
new production dependency added mid-session, not part of the pipeline's
original stage 6 scope.

**Why:** Without real coordinates, the override as designed in issue #23's
table ("gene falls inside a typed SCCmec cassette") literally cannot be
computed - there is no tool already in the pipeline that answers "where."
SCCmecExtractor is new (a 2026 preprint) but is pip/Docker packaged with
pinned releases (not just `:latest`) and reported 100% typing concordance
with `sccmec`/staphopia-sccmec-family tools on 1,454 *S. aureus* genomes in
its own benchmark - a reasonable, checked bet, not a blind one, but still an
unreviewed-in-a-journal tool worth revisiting if it turns out to be poorly
maintained.

**Affects:** Chromosomal tier (the SCCmec-override branch specifically).

**Author:** Claude (session with Daniel Garbozo)

## 2026-08-25: Tier-resolution rule engine gaps filled for stage 5f (issue #29)

**Decision:** The four-tier table in `CLAUDE.md`/issue #23 does not give an
explicit rule for every combination the five upstream signals (issues
#24-#28) can actually produce. `bin/tier_resolution.py`'s docstring has the
full precedence order and reasoning; the two real gaps and their resolution:

1. **Majority-chromosome classifier agreement** (2 of 3 classifiers say
   chromosome) has no defined tier - the table only names a "moderate"
   partial-confidence tier for majority-*plasmid* agreement, not the
   chromosomal direction. Resolved as: routes to Ambiguous, not Chromosomal.
2. **Plasmid classifier agreement (unanimous or majority) with zero
   mobility/replication markers** doesn't meet the floor either plasmid tier
   requires (both need ">=1 marker"), but isn't literally "classifiers
   disagree" either. Resolved as: routes to Ambiguous as the closest fit.

Precedence order (first match wins): SCCmec override (#28) beats everything,
including a length-floor or disagreement result that would otherwise apply -
direct structural evidence of cassette membership explains away classifier
confusion. Length floor (#27) is checked next, before classifier agreement
is even consulted, since a too-short contig makes that agreement untrustworthy
in the first place.

**Why:** Silently picking a default for an undefined table cell would be
exactly the kind of change this decisions.md file exists to prevent going
unnoticed. Both gaps were resolved toward the *more conservative* tier
(Ambiguous over a stronger claim) rather than guessing generously - this
pipeline's stated purpose is auditable, benchmarked confidence, not maximum
plasmid-call recall.

**Affects:** Ambiguous tier (both gaps), and the precedence relationship
between Chromosomal-via-override and every other tier. This is the actual
rule-engine logic stage 7 will validate - if real validation data shows
either gap-filling choice is wrong, that supersedes this entry, not a silent
code change.

**Author:** Claude (session with Daniel Garbozo)

## 2026-08-26: RFPlasmid channel join fix in evidence integration (no rule change)

**Decision:** Fixed a wiring bug in `subworkflows/local/evidence_integration/main.nf`:
`rfplasmid_prediction` was joined with a plain `.join()`, unlike `platon_tsv`
which already used `.join(remainder: true)`. Whenever RFPlasmid was skipped
(`--plasmid_skip_rfplasmid`, which the `test` profile sets), the plain join
produced an empty channel, so `CLASSIFIER_AGREEMENT` and everything
downstream of it - including `TIER_RESOLUTION` - silently received zero
inputs and never ran. Nextflow does not treat an empty channel as an error,
so the pipeline reported "completed successfully" with zero tier calls
produced. `bin/classifier_agreement.py`'s `--rfplasmid` argument was also
`required=True` (unlike `--platon`, already optional) and has been made
optional to match, since `agreement_level()` already computes agreement over
however many classifiers actually ran.

This is not a change to the tier rule table itself (CLAUDE.md's four-tier
table, or the precedence order from the 2026-08-25 tier-resolution entry) -
it is a bug that prevented the rule engine from running at all in any
configuration that skips RFPlasmid. Found by actually running
`-profile test,docker` end to end for the first time (prior local attempts
had failed earlier in the pipeline on unrelated issues - see git history for
the fARGene/SCCmecExtractor/CUSTOM_MULTIQC fixes from the same session).

**Why:** A rule engine that silently produces zero output on a valid,
documented pipeline configuration (the test profile explicitly sets
`plasmid_skip_platon` and `plasmid_skip_rfplasmid`, same as each other) is a
correctness bug, not a design choice - RFPlasmid was never meant to be
treated differently from Platon here.

**Affects:** No tier's rule definition changes. What changes is that
`TIER_RESOLUTION` now actually executes when 1 or 2 of the 3 classifiers are
unavailable, using whatever subset ran. Worth flagging separately: a
`unanimous_plasmid`/`unanimous_chromosome` call from `agreement_level()` does
not distinguish `n_classifiers == 1` from `n_classifiers == 3` - with only
MOB-suite running (as in `-profile test`), "unanimous" reflects one
classifier's vote, not the three-way agreement the tier table's
High-confidence/Moderate-confidence rows describe. `tier_resolution.py` does
not currently gate on `n_classifiers`, so results produced with Platon and/or
RFPlasmid skipped should not be reported as if all three classifiers agreed.
Not changed here since it is a rule-table design question, not a wiring bug -
flagging for a future decision once stage 7 (#22) is underway.

**Author:** Claude (session with Daniel Garbozo)
