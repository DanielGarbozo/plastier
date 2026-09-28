#!/usr/bin/env python3

"""
Stage 7b (issue #32): per-tier precision / recall / F1 against the closed-genome
ground truth.

Compares stage 5's `tier` call per ARG (bin/tier_resolution.py output) with where
that ARG actually sits in a closed reference genome.

Ground truth: assets/validation/closed_genome_ground_truth.csv lists, per
genome (`Genome_Accession`), every plasmid replicon's accession
(`Plasmid_Accession`). An ARG is truly on a plasmid iff the contig it was called
on is one of its genome's plasmid replicons, and truly chromosomal otherwise.
This is exact per-contig truth - matching on gene *names* instead would call a
gene like qacA "plasmid" everywhere in a genome that carries it on both a plasmid
and the chromosome.

That only works if `input_sequence_id` in the predictions is the reference
replicon's own accession, which is the case when stages 3-6 run directly on the
reference FASTA (see --assemblies) and NOT when reads are re-assembled (contigs
are then just "1", "2", ...). Predictions whose sequence ids are not accessions
are rejected rather than silently scored as all-chromosomal.

Metrics are computed per tier and never pooled across tiers, so a tier that is
systematically unreliable cannot hide behind a good aggregate (see docs/decisions.md).

- High / Moderate-confidence plasmid: positive class is "true plasmid".
- Chromosomal: positive class is "true chromosome".
- Ambiguous: a deliberate abstention with no true-class equivalent, so only the
  call counts are reported (precision/recall/F1 are left blank).

Per tier: TP = calls of that tier whose true location matches the tier's class,
FP = calls of that tier that do not, FN = ARGs of the tier's true class that were
NOT assigned that tier (so recall = share of all true-plasmid/true-chromosome
ARGs that ended up in this tier).
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

__version__ = "0.2.0"

TIER_SUFFIX = ".tier_resolution.tsv"
ACCESSION_RE = re.compile(r"^[A-Z]{1,2}_?[A-Z0-9]+\.\d+$")

PLASMID_TIERS = ("High-confidence plasmid", "Moderate-confidence plasmid")
CHROMOSOMAL_TIER = "Chromosomal"
AMBIGUOUS_TIER = "Ambiguous"
TIERS = (*PLASMID_TIERS, CHROMOSOMAL_TIER, AMBIGUOUS_TIER)


def first_token(seq_id) -> str:
    return str(seq_id).split()[0]


def genome_id_from_path(path: Path) -> str:
    name = path.name
    return name[: -len(TIER_SUFFIX)] if name.endswith(TIER_SUFFIX) else path.stem


def load_plasmid_replicons(path: str) -> dict:
    gt = pd.read_csv(path, dtype=str)
    missing = {"Genome_Accession", "Plasmid_Accession"} - set(gt.columns)
    if missing:
        sys.exit(f"Error: ground truth is missing column(s) {sorted(missing)}; found {list(gt.columns)}")
    gt = gt.dropna(subset=["Genome_Accession", "Plasmid_Accession"])
    return gt.groupby("Genome_Accession")["Plasmid_Accession"].agg(set).to_dict()


def label_predictions(pred_path: Path, genome_id: str, replicons: dict) -> pd.DataFrame:
    if genome_id not in replicons:
        sys.exit(f"Error: genome '{genome_id}' (from {pred_path.name}) is not in the ground truth; pass --genome-id.")

    df = pd.read_csv(pred_path, sep="\t", dtype=str)
    missing = {"gene_symbol", "input_sequence_id", "tier"} - set(df.columns)
    if missing:
        sys.exit(f"Error: {pred_path.name} is missing column(s) {sorted(missing)}.")
    if df.empty:
        return df.assign(genome_id=genome_id, true_location=[])

    df = df.copy()
    df["contig"] = df["input_sequence_id"].map(first_token)
    not_accessions = sorted(set(df["contig"][~df["contig"].str.match(ACCESSION_RE)]))
    if not_accessions:
        sys.exit(
            f"Error: {pred_path.name} has sequence ids that are not replicon accessions "
            f"({not_accessions[:5]}...). Per-contig truth needs stages 3-6 run on the reference "
            "FASTA (--assemblies), not on a re-assembly."
        )

    df["genome_id"] = genome_id
    df["true_location"] = df["contig"].map(lambda c: "Plasmid" if c in replicons[genome_id] else "Chromosome")
    return df


def tier_metrics(eval_df: pd.DataFrame) -> pd.DataFrame:
    n_true = eval_df["true_location"].value_counts()
    rows = []
    for tier in TIERS:
        calls = eval_df[eval_df["tier"] == tier]
        row = {
            "Tier": tier,
            "Total_Calls": len(calls),
            "Calls_True_Plasmid": int((calls["true_location"] == "Plasmid").sum()),
            "Calls_True_Chromosome": int((calls["true_location"] == "Chromosome").sum()),
            "TP": None, "FP": None, "FN": None, "Precision": None, "Recall": None, "F1_Score": None,
        }
        if tier != AMBIGUOUS_TIER:
            target = "Plasmid" if tier in PLASMID_TIERS else "Chromosome"
            tp = int((calls["true_location"] == target).sum())
            fp = len(calls) - tp
            fn = int(n_true.get(target, 0)) - tp
            precision = tp / (tp + fp) if tp + fp else 0.0
            recall = tp / (tp + fn) if tp + fn else 0.0
            f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
            row.update(TP=tp, FP=fp, FN=fn, Precision=round(precision, 4), Recall=round(recall, 4), F1_Score=round(f1, 4))
        rows.append(row)
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(
        prog="evaluate_metrics",
        description="Per-tier precision/recall/F1 of tier_resolution.py output against the closed-genome ground truth.",
    )
    parser.add_argument("--predictions", required=True, nargs="+", type=Path,
                        help="One or more *.tier_resolution.tsv files (one per genome)")
    parser.add_argument("--ground-truth", required=True, help="assets/validation/closed_genome_ground_truth.csv")
    parser.add_argument("--genome-id", help="Genome accession; only valid with a single predictions file. "
                                             f"Default: the filename minus '{TIER_SUFFIX}'.")
    parser.add_argument("--output", required=True, help="Output per-tier metrics TSV")
    parser.add_argument("--per-arg-output", help="Optional TSV of every ARG with its predicted tier and true location")
    parser.add_argument("--version", action="version", version=f"evaluate_metrics {__version__}")
    args = parser.parse_args()

    if args.genome_id and len(args.predictions) != 1:
        sys.exit("Error: --genome-id can only be used with a single predictions file.")

    replicons = load_plasmid_replicons(args.ground_truth)
    labelled = [
        label_predictions(p, args.genome_id or genome_id_from_path(p), replicons) for p in args.predictions
    ]
    eval_df = pd.concat(labelled, ignore_index=True)
    if eval_df.empty:
        sys.exit("Error: no ARG calls to evaluate.")

    tier_metrics(eval_df).to_csv(args.output, sep="\t", index=False)
    if args.per_arg_output:
        eval_df[["genome_id", "gene_symbol", "contig", "tier", "true_location"]].to_csv(
            args.per_arg_output, sep="\t", index=False
        )
    print(f"Evaluated {len(eval_df)} ARG calls across {eval_df['genome_id'].nunique()} genome(s) -> {args.output}")


if __name__ == "__main__":
    main()
