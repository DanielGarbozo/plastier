# Phase 3 scale-up: candidate accessions

About 20 public *S. aureus* short-read runs for the Phase 3 discovery run (see
[`roadmap.md`](roadmap.md)). The `--sra_ids` file is
[`assets/phase3/sra_ids_candidates.txt`](../assets/phase3/sra_ids_candidates.txt) (one run
per line, the format `workflows/plastier.nf` consumes); full ENA metadata is in
[`assets/phase3/candidates_metadata.tsv`](../assets/phase3/candidates_metadata.tsv).

**This does not set `sra_ids` in `conf/full.config`.** As with the pilot (issue #21), the list
needs a review before it is wired in. Until then it is passed explicitly with `--sra_ids`.

## How they were chosen

[`assets/phase3/select_candidates.py`](../assets/phase3/select_candidates.py) applies every
rule below with a fixed seed (42) against ENA's portal API; nothing was picked by hand.

1. ENA `scientific_name` is exactly *Staphylococcus aureus*; WGS, genomic, Illumina, paired-end.
2. 250-450 Mb of sequence (about 85-150x coverage) at 120-151 bp mean read length.
3. Country, collection date, host and isolation source are all recorded.
4. Not a Phase 0 pilot run, and not from a study a pilot run came from.
5. 15 human-host runs, at most one per study, spread over as many countries as possible,
   plus one run each from cat, cattle, dog, horse and pig. Companion and livestock hosts are
   included on purpose: plasmid-borne resistance is often host-associated.

ENA keeps growing, so re-running the script later can return different runs. This file and
`candidates_metadata.tsv` are the record of what was selected.

## Candidates (20 runs, 19 studies, 4.1 GB of FASTQ in total)

| # | Run | Host | Country | Collected | Isolation source | Sequence | FASTQ | Study |
| - | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | `SRR23309788` | Homo sapiens | USA: Pennsylvania | 2019-12 | throat | 313 Mb | 176 MB | PRJNA918392 |
| 2 | `SRR3731507` | Homo sapiens | United Kingdom: Brighton | 2011-12-13 | nose | 366 Mb | 177 MB | PRJNA324190 |
| 3 | `SRR12333138` | Homo sapiens | Egypt: Alexandria | 2015 | blood | 426 Mb | 252 MB | PRJNA648411 |
| 4 | `SRR2057031` | Homo sapiens | Australia: Northern Territory | 2013-10-03 | nasal | 434 Mb | 329 MB | PRJNA286158 |
| 5 | `DRR456378` | Homo sapiens | Japan:Gunma | 2020-03-09 | blood culture | 315 Mb | 183 MB | PRJDB15501 |
| 6 | `SRR27322154` | Homo sapiens | Spain: Sabadell, Barcelona | 2018 | blood | 282 Mb | 194 MB | PRJNA1055690 |
| 7 | `SRR21098763` | Homo sapiens | Denmark:Copenhagen | 2018 | blood | 256 Mb | 170 MB | PRJNA870263 |
| 8 | `SRR26374374` | Homo sapiens | Russia:Tomsk | 2021 | blood | 294 Mb | 122 MB | PRJNA1027597 |
| 9 | `SRR6332104` | Homo sapiens | New Zealand: Porirua | 2014 | Human clinical sample | 426 Mb | 254 MB | PRJNA412108 |
| 10 | `SRR6347285` | Homo sapiens | Cambodia: Siem Reap | 2008 | clinical sample | 444 Mb | 157 MB | PRJNA418899 |
| 11 | `SRR14862042` | Homo sapiens | Switzerland: St.Gallen | 2018-04-20 | skin | 374 Mb | 289 MB | PRJNA738326 |
| 12 | `SRR7867492` | Homo sapiens | China: Shaanxi | 2011-06-08 | Hospital patients | 442 Mb | 218 MB | PRJNA491594 |
| 13 | `SRR36881360` | Homo sapiens | South Korea: Chungcheong-buk-do, Cheongju-si | 2001 | wound | 340 Mb | 234 MB | PRJNA1387400 |
| 14 | `SRR18339536` | Homo sapiens | Mexico: Mexico City | 2018-04-04 | Right-Popliteal-Fossa | 417 Mb | 256 MB | PRJNA816913 |
| 15 | `SRR14267647` | Homo sapiens | France | 2015 | sputum from cystic fibrosis patient | 279 Mb | 171 MB | PRJNA721116 |
| 16 | `DRR484283` | Felis catus | Japan:Hokkaido | 2020-03-20 | catheteruria | 291 Mb | 159 MB | PRJDB16052 |
| 17 | `SRR29758699` | Bos taurus | New Zealand: North Canterbury | 2021-11-30 | Bacterial Isolate from bulk tank milk sample | 356 Mb | 196 MB | PRJNA1130542 |
| 18 | `DRR484144` | Canis lupus familiaris | Japan:Hiroshima | 2018-11-06 | scrotal mass | 263 Mb | 142 MB | PRJDB16052 |
| 19 | `SRR13450679` | Equus caballus | Switzerland | 2019-07 | pastern skin | 319 Mb | 151 MB | PRJNA692738 |
| 20 | `ERR7434454` | Sus scrofa domesticus | Spain | 2013 | Carrier pig | 430 Mb | 282 MB | PRJEB49018 |

## What this list does and does not establish

- **The species label is ENA's, not independently verified.** The pilot list once turned out to
  be 7/10 wrong. Here `scientific_name` is filtered exactly and the runs are ordinary
  isolate WGS, but a mislabelled or mixed sample would still get through. The pipeline gives
  an independent check for free: stage 6's MLST should return a *S. aureus* sequence type
  for every sample. A sample with no ST should be treated as a bad accession, not a result.
- **No strain typing is claimed.** ST, spa and SCCmec type are not known until stage 6 runs.
- Four countries appear more than once (New Zealand, Japan, Spain, Switzerland), because the human and animal picks are made independently; the 15 human runs are from 15 different countries. The cat and dog runs share one study (PRJDB16052): each animal host has only one to six qualifying studies in ENA, so the one-run-per-study rule holds for the human runs but not across hosts.
- Sequencing depth is 85-150x by design; results should not be read as an effect of depth
  (that is what the stage 7 sensitivity analysis, #35, is for).
