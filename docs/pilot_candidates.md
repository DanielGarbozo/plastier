# Phase 0 pilot: candidate accessions (issue #21)

Candidate shortlist for `conf/pilot.config`'s `sra_ids`, per [issue #21](https://github.com/DanielGarbozo/plastier/issues/21).
The matching `--sra_ids` file is [`assets/pilot/sra_ids_candidates.csv`](../assets/pilot/sra_ids_candidates.csv)
(one run accession per line, per the format `workflows/plastier.nf` actually consumes).

**This does not edit `conf/pilot.config`.** Per the issue: "please don't edit
`conf/pilot.config` directly - the final list needs review before it's wired in."

## Important: an existing list in `conf/pilot.config` needs to be reverted

Commit `8ee1803` ("Agrega sra_ids a pilot.config", Rodrigo Puertas, 2026-08-25) set
`sra_ids` directly on `conf/pilot.config`, bypassing the review process this issue
asks for. Verifying those 10 accessions against ENA's public API (`instrument_platform`,
`scientific_name` fields) found that **7 of the 10 are not *S. aureus*** at all:

| Accession | `scientific_name` (from ENA) |
| --- | --- |
| ERR017635 | *Staphylococcus aureus* ✅ |
| ERR018278 | *Sarcophilus harrisii* (Tasmanian devil) ❌ |
| ERR018282 | *Homo sapiens*, single-end ❌ |
| SRR1635435 | *Homo sapiens*, single-end ❌ |
| ERR103394 | *Staphylococcus aureus* ✅ |
| SRR5244198 | *Homo sapiens*, single-end ❌ |
| ERR340820 | *Streptococcus pneumoniae* ❌ |
| SRR7973517 | *Mus musculus* ❌ |
| ERR410049 | *Staphylococcus aureus* ✅ |
| ERR245594 | *Danio rerio* ❌ |

**Recommendation:** revert `sra_ids` on `conf/pilot.config` back to `null` (its
original, intentional fail-closed state) until a reviewed list - such as the one
below - is merged. The 3 valid *S. aureus* accessions (`ERR017635`, `ERR103394`,
`ERR410049`) are real but have no confirmed strain/ST/spa/SCCmec metadata, so they
are not included in the candidate list below; they could be investigated further
in a follow-up if useful.

## Candidate list (9 of 10 slots filled)

All below verified against ENA's public API (`instrument_platform=ILLUMINA`,
`library_layout=PAIRED`, real `scientific_name=Staphylococcus aureus`). File sizes
are the smaller of multiple public runs per strain where more than one existed, to
keep pilot-scale compute modest per the issue's request.

| # | Accession | Strain | MRSA/MSSA | ST | spa | SCCmec | Approx. size (R1+R2) | Source |
| - | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `SRR13251326` | MW2 | MRSA (community-acquired) | ST1 | - | IVa | ~97 MB | Baba et al. 2002, reference genome; ENA study PRJNA685119 |
| 2 | `SRR13251228` | USA300 LAC | MRSA (community-acquired) | ST8 | - | IVa | ~100 MB | Diep et al., USA300 lineage; ENA study PRJNA685119 |
| 3 | `SRR13251312` | MRSA252 | MRSA (hospital-acquired, EMRSA-16) | ST36 | - | II | ~124 MB | Holden et al. 2004 (PNAS); ENA study PRJNA685119 |
| 4 | `SRR13251283` | N315 | MRSA (hospital-acquired, historic) | ST5 | - | II | ~131 MB | Kuroda et al. 2001, first sequenced MRSA genome; ENA study PRJNA685119 |
| 5 | `ERR046448` | TW20 | MRSA (hospital outbreak) | ST239 | - | - | ~92 MB | Sanger Institute / Holden et al., UK ICU outbreak strain; ENA study PRJEB2627 |
| 6 | `SRR11782467` | MSSA476 | MSSA | ST1 | - | N/A | ~48 MB | Holden et al. 2004 (PNAS), paired with MRSA252 in same study; ENA study PRJNA632444 |
| 7 | `SRR1948056` | NCTC 8325 | MSSA | ST8 | - | N/A | ~63 MB | Gillaspy et al., classic reference genome; ENA study PRJNA279815 |
| 8 | `ERR16669663` | Newman | MSSA | ST8 | - | N/A | ~95 MB | Baba et al. 2008, classic clinical isolate; ENA study PRJEB72799 |
| 9 | `SRR9141411` | GD487 | MSSA | ST398 | - | N/A | ~140 MB | High-virulence, human-associated ST398 MSSA; ENA study PRJNA541965 |
| 10 | *(open)* | - | - | - | - | - | - | See note below |

**Diversity achieved:** 5 distinct STs on the MRSA side (1, 8, 36, 5, 239), 3 distinct
STs on the MSSA side (1, 8, 398) - though NCTC 8325 and Newman share ST8 (different
isolates, different history, but same MLST background; worth knowing when
interpreting strain-typing-stage results). A genuine mix of MRSA/MSSA and of
SCCmec types (II, IVa) is present.

## Slot 10: intentionally left open

A 10th candidate (originally RN4220, a common lab strain) was investigated and
dropped: every public Illumina paired-end run found for RN4220 carries an
artificially introduced plasmid or transposon construct for unrelated molecular
biology experiments (e.g. `RN4220/pG0400`, `RN4220::Tn5405` - see ENA studies
PRJNA289526 and PRJEB2655). Using one of these would make `plastier` call "has a
plasmid" on a cloning vector rather than real biology, which would corrupt the
tier-resolution validation this pilot exists to support. Rather than force a
weaker substitute under time pressure, this slot is left for the team to fill
after review - happy to keep looking if useful.

## Caveats

- Approximate sizes above are `fastq_bytes` (R1+R2) as reported by ENA at time of
  writing; actual download size may vary slightly.
- `spa` type is left blank throughout except where a specific PubMLST/paper source
  confirmed it (only HO 5096 0412 had a confirmed `spa` type in the sources found,
  and that strain was ultimately replaced by N315 for lack of a public Illumina
  paired-end run - see above).
- These are reference/type strains long used in *S. aureus* genomics, not a
  clinical surveillance cohort - appropriate for pipeline validation (the ST/SCCmec
  ground truth is very well established in the literature), but not representative
  of current circulating epidemiology.
