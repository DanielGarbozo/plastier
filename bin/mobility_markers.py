#!/usr/bin/env python3

"""
Stage 5b (issue #25): replicon / mobility marker evidence.

Reports, per ARG, whether its contig carries a replicon-typing or
relaxase/mobility marker per MOB-suite's own contig_report.txt - the second
required signal for the "High-confidence plasmid" and "Moderate-confidence
plasmid" tiers in the four-tier framework (CLAUDE.md, issue #23): both
require "≥1 mobility/replication marker" on top of classifier agreement.

Sourced from MOB-suite alone (no new tool) - contig_report.txt already
carries `rep_type(s)` (replicon typing) and `relaxase_type(s)` (conjugative
relaxase typing) per contig, plus MOB-suite's own `predicted_mobility` call
(conjugative/mobilizable/non-mobilizable/-). MOB-suite prints '-' for "no
hit", not a blank cell, and can report multiple types separated by commas -
both are treated as "no marker" / "has a marker" respectively without
needing to parse the individual type list.

Same join-key and sample-scoping approach as bin/classifier_agreement.py
(issue #24, see that file's docstring for why): join on the first
whitespace-delimited token of the contig ID, and filter the merged,
all-samples hAMRonization report to --sample-id's own rows.
"""

import argparse

import pandas as pd

__version__ = "0.1.0"

MARKER_COLUMNS = ("rep_type(s)", "relaxase_type(s)")


def first_token(contig_id: str) -> str:
    return str(contig_id).split()[0]


def is_present(value) -> bool:
    value = str(value).strip()
    return bool(value) and value != "-"


def load_hamronization(path: str, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["input_file_name"] == sample_id].copy()
    df["contig_key"] = df["input_sequence_id"].map(first_token)
    return df


def load_mobsuite_markers(path: str) -> dict:
    """contig_key -> {rep_type, relaxase_type, predicted_mobility, has_marker}"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["contig_key"] = df["contig_id"].map(first_token)
    df["has_marker"] = df[list(MARKER_COLUMNS)].apply(lambda row: any(is_present(v) for v in row), axis=1)
    df = df.rename(columns={"rep_type(s)": "rep_type", "relaxase_type(s)": "relaxase_type"})
    return df.set_index("contig_key")[["rep_type", "relaxase_type", "predicted_mobility", "has_marker"]].to_dict("index")


def main():
    parser = argparse.ArgumentParser(
        prog="mobility_markers",
        description="Stage 5b: per-ARG replicon/mobility marker evidence from MOB-suite.",
    )
    parser.add_argument("--hamronization", required=True, help="hAMRonization merged ARG report (tsv, all samples)")
    parser.add_argument("--sample-id", required=True, help="input_file_name value identifying this sample's rows")
    parser.add_argument("--mobsuite", required=True, help="MOB-suite contig_report.txt")
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"mobility_markers {__version__}")
    args = parser.parse_args()

    hamronization = load_hamronization(args.hamronization, args.sample_id)
    markers = load_mobsuite_markers(args.mobsuite)

    rows = []
    for _, arg_row in hamronization.iterrows():
        contig_key = arg_row["contig_key"]
        marker = markers.get(contig_key, {})
        rows.append(
            {
                "gene_symbol": arg_row.get("gene_symbol"),
                "input_sequence_id": arg_row["input_sequence_id"],
                "contig_key": contig_key,
                "rep_type": marker.get("rep_type"),
                "relaxase_type": marker.get("relaxase_type"),
                "predicted_mobility": marker.get("predicted_mobility"),
                "has_marker": marker.get("has_marker", False),
            }
        )

    columns = [
        "gene_symbol",
        "input_sequence_id",
        "contig_key",
        "rep_type",
        "relaxase_type",
        "predicted_mobility",
        "has_marker",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
