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
  plasmid_(X).fasta file. Rows are kept only where molecule_type is "plasmid"
  (chromosome-assigned contigs are excluded, matching what a plasmid bin means).

  Three details that only show up on a real MOB-suite report (found running the
  pipeline on the closed genomes; a hand-built fixture hid all three):
  - molecule_type is lowercase ("plasmid"), so the match is case-insensitive.
  - contig_id is the whole FASTA header ("NZ_LR027879.1 Staphylococcus aureus ...");
    PlasEval matches contigs by exact id, so only its first token is used - the
    replicon accession, which is what the ground truth uses.
  - primary_cluster_id is a MOB-cluster name (e.g. AA849) shared by unrelated
    plasmids in different genomes, and is "-" for a plasmid contig with no cluster.
    Bins are therefore named <sample_id>|<cluster>, and an unclustered contig is a
    bin of its own (<sample_id>|<contig>) instead of all sharing one "-" bin.
  Several reports (one per genome) can be converted into one file, since contig
  accessions are unique across genomes.

Usage:
    python plaseval_convert.py ground-truth \\
        --input assets/validation/closed_genome_ground_truth.csv \\
        --output gt_bins.tsv

    python plaseval_convert.py predicted \\
        --input results/mobsuite/*/contig_report.txt \\
        --output pred_bins.tsv
"""

import argparse

import pandas as pd

__version__ = "0.2.0"


def convert_ground_truth(input_path: str, output_path: str, only: list = None) -> None:
    gt = pd.read_csv(input_path, dtype=str)
    if only:  # a partial run must not be penalised for genomes it never processed
        gt = gt[gt["Genome_Accession"].isin(only)]
    plasmids = gt[gt["Plasmid_Accession"] != "None"].copy()
    out = pd.DataFrame(
        {
            "plasmid": plasmids["Plasmid_Accession"],
            "contig": plasmids["Plasmid_Accession"],
            "contig_len": plasmids["Plasmid_Size_bp"],
        }
    )
    out.to_csv(output_path, sep="\t", index=False)


def convert_predicted(input_paths: list, output_path: str) -> None:
    frames = []
    for path in input_paths:
        mob = pd.read_csv(path, sep="\t", dtype=str)
        plasmids = mob[mob["molecule_type"].str.strip().str.lower() == "plasmid"]
        contigs = plasmids["contig_id"].str.split().str[0]
        cluster = plasmids["primary_cluster_id"].fillna("-").str.strip()
        unclustered = cluster.isin(["", "-"])
        bin_name = plasmids["sample_id"] + "|" + cluster.where(~unclustered, contigs)
        frames.append(pd.DataFrame({"plasmid": bin_name, "contig": contigs, "contig_len": plasmids["size"]}))
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=["plasmid", "contig", "contig_len"])
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
    gt_parser.add_argument("--only", nargs="+", help="Keep only these Genome_Accession values (partial runs)")

    pred_parser = sub.add_parser("predicted", help="Convert MOB-suite's contig_report.txt to pred_bins.tsv")
    pred_parser.add_argument("--input", required=True, nargs="+",
                             help="One or more contig_report.txt files (stage 4, MOB-suite), one per genome")
    pred_parser.add_argument("--output", required=True, help="Output path for pred_bins.tsv")

    parser.add_argument("--version", action="version", version=f"plaseval_convert {__version__}")
    args = parser.parse_args()

    if args.mode == "ground-truth":
        convert_ground_truth(args.input, args.output, args.only)
    elif args.mode == "predicted":
        convert_predicted(args.input, args.output)


if __name__ == "__main__":
    main()
