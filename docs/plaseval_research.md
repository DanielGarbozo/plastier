# PlasEval research and integration (issue #33)

Per [issue #33](https://github.com/DanielGarbozo/plastier/issues/33)'s explicit
first step: *"Research PlasEval's actual input format/requirements (not yet
done in this repo - start here before writing any integration code)."* This
covers that research, plus a tested converter and a documented methodological
gap that needs a decision before this can run on real pipeline output.

## What PlasEval is

[PlasEval](https://github.com/cchauve/PlasEval) (Mane, Sanderson, White,
Zaheer, Beiko, Chauve; *BMC Bioinformatics* 2024,
[10.1186/s12859-024-05941-0](https://doi.org/10.1186/s12859-024-05941-0))
evaluates **plasmid binning** - whether contigs from a draft assembly were
correctly grouped into the individual plasmids they actually belong to. This
is a different question from #32's per-tier precision/recall, which only asks
"is this ARG's contig a plasmid or the chromosome" - PlasEval asks "if there
are two separate plasmids in this genome, did the pipeline correctly tell
their contigs apart, or did it merge/split them incorrectly." Confirmed by
cloning the tool and running it directly (see "Tested" below), not just read
from its docs.

PlasEval has two modes: `eval` (predicted bins vs. a ground truth - precision/
recall/F1) and `comp` (dissimilarity between any two bin sets, e.g. two
different tools). Issue #33 only asks for `eval` mode.

## Input format

Both modes take plasmid-bins files: a TSV with a header row and exactly three
columns:

```
plasmid	contig	contig_len
```

One row per (bin, contig) pair - `plasmid` is the bin's identifier, `contig`
the contig's identifier, `contig_len` its length in bp. Confirmed directly
from PlasEval's own README and its bundled `examples/input/*.tsv` files.

Command:

```sh
python plaseval.py eval --pred PREDICTED_BINS_TSV --gt GROUNDTRUTH_BINS_TSV \
    --out_file OUT_FILE --log_file LOG_FILE (--min_len LEN_THRESHOLD)
```

## Converter: `assets/validation/plaseval_convert.py`

Converts plastier's own data into that three-column format. Two modes:

- **`ground-truth`**: `assets/validation/closed_genome_ground_truth.csv` (#31)
  -> `gt_bins.tsv`. Each ground-truth plasmid row already represents one
  whole closed replicon, so `plasmid == contig == Plasmid_Accession` and
  `contig_len == Plasmid_Size_bp`.
- **`predicted`**: a stage-4 MOB-suite `contig_report.txt` -> `pred_bins.tsv`.
  MOB-recon's own `primary_cluster_id` column groups contigs into the same
  reconstructed plasmid (confirmed against mob-suite's source,
  `mob_suite/constants.py`'s `MOB_RECON_INFO_HEADER` - the README's own
  documentation table is missing the `contig_id` column that the actual code
  writes, so the source was checked directly rather than trusting the docs
  alone). Rows are kept only where `molecule_type == "Plasmid"`.

### Tested

- `ground-truth` mode run against the real, already-merged
  `closed_genome_ground_truth.csv`: converts cleanly, 42 plasmid rows in, 42
  out (1 genome in that set, Mu50-style chromosome-only genomes, if any get
  added later - see `docs/stage7_ground_truth.md` - correctly contribute zero
  rows).
- `predicted` mode run against a hand-built, clearly-synthetic
  `contig_report.txt` matching MOB-suite's real confirmed column schema:
  correctly filters out the chromosome row and groups by
  `primary_cluster_id`.
- The actual `plaseval.py eval` command run end-to-end (PlasEval cloned and
  installed in this environment, not just read about) against converter
  output: when contig IDs are shared between predicted and ground truth, it
  returns exactly the expected result (precision = recall = F1 = 1.0 for a
  perfect match). This confirms the converter's output is valid PlasEval
  input, not just correctly-shaped.

## Open problem: contig identity (read before running this for real)

PlasEval's precision/recall is computed from **shared contig identifiers**
between the predicted and ground-truth bin files - it does not fuzzy-match by
length or content. Confirmed directly: converting a synthetic
`contig_report.txt` whose `contig_id` values were assembler-style names
(`contig_2`, `contig_3`) against the real ground truth (whose contig IDs are
NCBI accessions, e.g. `NZ_CP149493.1`) scored **0.0 precision and recall
across the board**, even though the contig lengths matched exactly. Re-running
the identical comparison with the contig IDs changed to match the ground
truth's accessions scored a perfect 1.0.

This means: if stage 4 runs on genomes assembled from scratch (fresh
`contig_1`, `contig_2`, ... names from Unicycler), PlasEval will show zero
overlap with the ground truth even when the pipeline's actual plasmid content
is correct, because the contig names themselves never match. **The same issue
was designed around for #32's `tier_benchmark.py`** by running stages 3-6
directly on each ground-truth genome's own reference FASTA (skipping
assembly), so `input_sequence_id` is the original NCBI accession rather than
a freshly generated contig name. The same approach should be used here for
consistency between #32 and #33's results - running MOB-suite directly on the
reference replicons, not on a fresh assembly of simulated/real reads.

## Not done here (needs real pipeline output)

Actually running this against plastier's own stage 4 output requires running
MOB-suite (or the full pipeline) on the #31 ground-truth genomes - which
needs compute this environment doesn't have. Once that output exists
(ideally produced the way described above, replicon accessions preserved as
contig IDs), running the benchmark is two conversions plus one command:

```sh
python assets/validation/plaseval_convert.py ground-truth \
    --input assets/validation/closed_genome_ground_truth.csv --output gt_bins.tsv
python assets/validation/plaseval_convert.py predicted \
    --input <path/to/contig_report.txt> --output pred_bins.tsv
python plaseval.py eval --pred pred_bins.tsv --gt gt_bins.tsv \
    --out_file plaseval_report.tsv --log_file plaseval.log
```
