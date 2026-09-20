#!/usr/bin/env python3

"""
Stage 7c (issue #33): converts plastier's own data into the plasmid-bins TSV
format PlasEval (https://github.com/cchauve/PlasEval) expects for its `eval`
mode - three columns, tab-separated, with a header row: plasmid, contig,
contig_len. See docs/plaseval_research.md for how this format and the
`primary_cluster_id` mapping below were confirmed.

Two conversions:

  ground-truth   closed_genome_ground_truth.csv (#31)  -> gt_bins.tsv
  predicted      MOB-suite's contig_report.txt (stage 4) -> pred_bins.tsv

Ground truth mapping:
  Each non-chromosome row in the ground-truth CSV already represents one
  whole, closed plasmid replicon - there is no sub-contig structure to
  reconstruct. [Judgement call] So each ground-truth "bin" is exactly one
  replicon: plasmid == contig == Plasmid_Accession, contig_len ==
  Plasmid_Size_bp. Rows with Plasmid_Accession == "None" (chromosome-only
  genomes, e.g. Mu50 if added - see docs/stage7_ground_truth.md) contribute
  no rows here, since PlasEval bins describe plasmid content only.

Predicted-bins mapping (from MOB-suite's contig_report.txt):
  [Judgement call] MOB-recon's own clustering identifier, primary_cluster_id
  ("primary MOB-cluster id of neighbor" per mob-suite's own README/
  constants.py), is what groups separate contigs into the same reconstructed
  plasmid - the same grouping MOB-recon itself uses to write each
  plasmid_(X).fasta file. That is used directly as PlasEval's `plasmid` bin
  ID. Rows are kept only where molecule_type == "Plasmid" (chromosome-
  assigned contigs are excluded, matching what a plasmid bin means).

Usage:
    python plaseval_convert.py ground-truth \\
        --input assets/validation/closed_genome_ground_truth.csv \\
        --output gt_bins.tsv

    python plaseval_convert.py predicted \\
        --input results/mobsuite/SAMPLE_contig_report.txt \\
        --output pred_bins.tsv
"""

import argparse

import pandas as pd

__version__ = "0.1.0"


def convert_ground_truth(input_path: str, output_path: str) -> None:
    gt = pd.read_csv(input_path, dtype=str)
    plasmids = gt[gt["Plasmid_Accession"] != "None"].copy()
    out = pd.DataFrame(
        {
            "plasmid": plasmids["Plasmid_Accession"],
            "contig": plasmids["Plasmid_Accession"],
            "contig_len": plasmids["Plasmid_Size_bp"],
        }
    )
    out.to_csv(output_path, sep="\t", index=False)


def convert_predicted(input_path: str, output_path: str) -> None:
    mob = pd.read_csv(input_path, sep="\t", dtype=str)
    plasmid_contigs = mob[mob["molecule_type"] == "Plasmid"].copy()
    out = pd.DataFrame(
        {
            "plasmid": plasmid_contigs["primary_cluster_id"],
            "contig": plasmid_contigs["contig_id"],
            "contig_len": plasmid_contigs["size"],
        }
    )
    out.to_csv(output_path, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        prog="plaseval_convert",
        description="Stage 7c: convert plastier data to PlasEval's plasmid-bins TSV format.",
    )
    sub = parser.add_subparsers(dest="mode", required=True)

    gt_parser = sub.add_parser("ground-truth", help="Convert closed_genome_ground_truth.csv (#31) to gt_bins.tsv")
    gt_parser.add_argument("--input", required=True, help="Path to closed_genome_ground_truth.csv")
    gt_parser.add_argument("--output", required=True, help="Output path for gt_bins.tsv")

    pred_parser = sub.add_parser("predicted", help="Convert MOB-suite's contig_report.txt to pred_bins.tsv")
    pred_parser.add_argument("--input", required=True, help="Path to a contig_report.txt (stage 4, MOB-suite)")
    pred_parser.add_argument("--output", required=True, help="Output path for pred_bins.tsv")

    parser.add_argument("--version", action="version", version=f"plaseval_convert {__version__}")
    args = parser.parse_args()

    if args.mode == "ground-truth":
        convert_ground_truth(args.input, args.output)
    elif args.mode == "predicted":
        convert_predicted(args.input, args.output)


if __name__ == "__main__":
    main()
