# Stage 7 ground truth: closed-genome benchmark set (issue #22)

Per [issue #22](https://github.com/DanielGarbozo/plastier/issues/22)'s scope: *"Curate/acquire
a closed-genome ground-truth set for S. aureus (complete assemblies with resolved
plasmid content) - a separate concern from the Phase 0 pilot's raw short-read input
(issue #21)."*

This is a **first tranche of 3 genomes**, not the final set - see "What's not
covered yet" below. Machine-readable version:
[`assets/ground_truth/closed_genomes_stage7.tsv`](../assets/ground_truth/closed_genomes_stage7.tsv).

## Why these 3

All three are complete (closed) RefSeq/GenBank assemblies where every replicon
(chromosome vs. each plasmid) is a separately accessioned, finished sequence - so
"is this AMR gene on a plasmid or the chromosome" has a definitive, citable answer,
which is exactly what stage 5's tier calls need to be checked against.

They were chosen to cover three distinct ground-truth *shapes*, on purpose:

| Strain | Shape | Why it matters for validation |
| --- | --- | --- |
| **USA300 FPR3757** | Chromosomal AMR (*mecA*/SCCmec) **and** plasmid-borne AMR (pUSA03) in the same genome | Tests whether the pipeline can correctly split AMR calls between chromosome and plasmid *within one genome*, not just across genomes |
| **N315** | Chromosomal AMR (*mecA*/SCCmec) **and** plasmid-borne AMR (pN315, *blaZ*) | A second, independent instance of the "both" case, different SCCmec type and different plasmid/gene pair, from the strain that later became the reference against which the field's *S. aureus* nomenclature was built |
| **Mu50** | Chromosomal AMR only - **no** resistance plasmid | Negative control: any plasmid-tier call the pipeline makes on Mu50 is by definition a false positive, since the source literature explicitly confirms this strain carries no penicillinase/resistance plasmid (unlike its close relative N315) |

## Replicon-level detail

### USA300 FPR3757
- Chromosome: `CP000255` - carries *mecA* inside a SCCmec IV element.
- Plasmid pUSA01: `CP000256` (~3.1 kb).
- Plasmid pUSA02: `CP000257`.
- Plasmid pUSA03: `CP000258` (~27 kb) - carries a cadmium-resistance/antimicrobial-resistance operon. This is the plasmid-borne AMR signal for this genome.
- Source: Diep BA et al. 2006, *Lancet* 367(9512):731-9, "Complete genome sequence of USA300, an epidemic clone of community-acquired meticillin-resistant Staphylococcus aureus."

### N315
- Chromosome: `BA000018` - carries *mecA* inside a SCCmec II element.
- Plasmid pN315: `AP003139` - carries *blaZ* (penicillinase/beta-lactamase). Confirmed directly in NCBI record `NG_047533.1` ("N315 pN315 blaZ gene for penicillinase").
- Source: Kuroda M et al. 2001, *Lancet* 357(9264):1225-40, "Whole genome sequencing of meticillin-resistant Staphylococcus aureus."

### Mu50
- Chromosome: `BA000017.4` - carries *mecA* inside a SCCmec II element (same type as N315; Mu50 and N315 are closely related, both ST5).
- No plasmid. Explicitly confirmed in Katayama Y et al. 2016, *Antimicrob Agents Chemother*, "Complete Reconstitution of the Vancomycin-Intermediate Staphylococcus aureus Phenotype of Strain Mu50 in Vancomycin-Susceptible S. aureus": the paper describes constructing a *N315* derivative "in which ... the plasmid carrying the gene for penicillinase (PCase; beta-lactamase) was eliminated" specifically "to mimic Mu3 and Mu50, which have no PCase plasmid." Mu50's intermediate vancomycin resistance (VISA phenotype) comes from chromosomal mutations, not a plasmid.
- Source: Kuroda M et al. 2001 (same paper as N315 - both were sequenced and reported together).

## What's not covered yet

This is a starting set, chosen for citability and to cover the three most basic
ground-truth shapes (chromosome-only, plasmid-only-signal, and both-in-one-genome).
It does **not** yet cover:

- **MSSA genomes** (no resistance plasmid *and* no SCCmec - a second kind of true
  negative, distinct from Mu50's "SCCmec present, plasmid absent" case).
- **Multi-plasmid AMR genomes** where different plasmids carry *different* AMR
  genes (tests whether per-plasmid attribution, not just plasmid-vs-chromosome,
  works).
- Strains outside the ST5/ST8 lineages already used for the issue #21 pilot
  candidates - some deliberate non-overlap with the pilot's short-read strains
  would strengthen the benchmark by testing on genomes the pipeline hasn't
  effectively "seen" via a near-identical isolate already.
- Any genome where a resistance gene sits on an SCC*mec* cassette specifically
  (as opposed to elsewhere on the chromosome) - relevant to stage 5's SCCmec
  override signal (`sccmecoverride/`) and worth a dedicated case once the SCCmec
  boundary coordinates for one of these strains are confirmed.

Recommend treating this as an incremental PR - happy to keep expanding once this
first set is reviewed, rather than delay review while chasing a complete set.
