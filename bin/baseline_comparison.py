#!/usr/bin/env python3

"""
Stage 7d (issue #34): single-tool baseline vs the four-tier framework.

Most plasmid-classification pipelines report one tool's plasmid/chromosome call as
fact. This scores exactly that - MOB-suite alone, Platon alone, RFPlasmid alone -
against the same per-contig ground truth as bin/evaluate_metrics.py (#32), next to
the tier framework, so the benchmark can show whether combining classifiers and
evidence actually beats picking one.

Every method is reduced to a plasmid / chromosome call per ARG:
- Tier framework: High- or Moderate-confidence plasmid -> plasmid, Chromosomal ->
  chromosome, Ambiguous -> no call (it abstains by design).
- Single classifier: its own `*_call` column; no call if the classifier had none
  (e.g. it was skipped, or never saw the contig).

Abstentions are what make a naive comparison misleading: a method that only answers
when it is sure looks more precise. So each row also reports Coverage, the share of
all ARGs the method answered, and Recall is always over ALL true plasmid (or
chromosome) ARGs, so an abstention costs recall rather than disappearing.

Per method and class (Plasmid / Chromosome):
    TP = called this class and truly is,  FP = called this class and truly is not,
    FN = truly this class but not called it (wrong call OR no call).

These are the same TP/FP/FN definitions as evaluate_metrics.py, and numbers here
answer a different question from #32's per-tier table (which never pools tiers).
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from evaluate_metrics import CHROMOSOMAL_TIER, PLASMID_TIERS, genome_id_from_path, label_predictions, load_plasmid_replicons  # noqa: E402

__version__ = "0.1.0"

CLASSIFIER_COLUMNS = {
    "MOB-suite alone": "mobsuite_call",
    "Platon alone": "platon_call",
    "RFPlasmid alone": "rfplasmid_call",
}
CLASSES = ("Plasmid", "Chromosome")


def tier_call(tier: str):
    if tier in PLASMID_TIERS:
        return "Plasmid"
    if tier == CHROMOSOMAL_TIER:
        return "Chromosome"
    return None  # Ambiguous: abstains


def classifier_call(value):
    if pd.isna(value):
        return None
    return {"plasmid": "Plasmid", "chromosome": "Chromosome"}.get(str(value).strip().lower())


def method_calls(eval_df: pd.DataFrame) -> dict:
    """method name -> Series of 'Plasmid' / 'Chromosome' / None, aligned to eval_df."""
    methods = {"Tier framework (High+Moderate = plasmid)": eval_df["tier"].map(tier_call)}
    for name, column in CLASSIFIER_COLUMNS.items():
        # A classifier that was skipped has no data at all - leave it out rather than
        # report an all-abstain row that looks like a result.
        if column in eval_df.columns and eval_df[column].notna().any():
            methods[name] = eval_df[column].map(classifier_call)
    return methods


def score(calls: pd.Series, truth: pd.Series) -> list:
    rows = []
    coverage = round(calls.notna().mean(), 4)
    for cls in CLASSES:
        tp = int(((calls == cls) & (truth == cls)).sum())
        fp = int(((calls == cls) & (truth != cls)).sum())
        fn = int(((truth == cls) & (calls != cls)).sum())
        precision = tp / (tp + fp) if tp + fp else 0.0
        recall = tp / (tp + fn) if tp + fn else 0.0
        f1 = 2 * precision * recall / (precision + recall) if precision + recall else 0.0
        rows.append({
            "Class": cls, "Calls": int((calls == cls).sum()), "Coverage": coverage,
            "TP": tp, "FP": fp, "FN": fn,
            "Precision": round(precision, 4), "Recall": round(recall, 4), "F1_Score": round(f1, 4),
        })
    return rows


def main():
    parser = argparse.ArgumentParser(
        prog="baseline_comparison",
        description="Score single classifiers against the four-tier framework on the closed-genome ground truth.",
    )
    parser.add_argument("--predictions", required=True, nargs="+", type=Path,
                        help="One or more *.tier_resolution.tsv files (one per genome)")
    parser.add_argument("--ground-truth", required=True, help="assets/validation/closed_genome_ground_truth.csv")
    parser.add_argument("--output", required=True, help="Output comparison TSV")
    parser.add_argument("--version", action="version", version=f"baseline_comparison {__version__}")
    args = parser.parse_args()

    replicons = load_plasmid_replicons(args.ground_truth)
    eval_df = pd.concat(
        [label_predictions(p, genome_id_from_path(p), replicons) for p in args.predictions], ignore_index=True
    )
    if eval_df.empty:
        sys.exit("Error: no ARG calls to evaluate.")

    rows = []
    for method, calls in method_calls(eval_df).items():
        for row in score(calls, eval_df["true_location"]):
            rows.append({"Method": method, **row})
    pd.DataFrame(rows).to_csv(args.output, sep="\t", index=False)
    print(f"Compared {len(rows) // len(CLASSES)} method(s) on {len(eval_df)} ARG calls "
          f"across {eval_df['genome_id'].nunique()} genome(s) -> {args.output}")


if __name__ == "__main__":
    main()
