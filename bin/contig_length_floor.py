#!/usr/bin/env python3

"""
Stage 5d (issue #27): validated contig-length floor.

Reports, per ARG, whether its contig's length falls below a validated
minimum, below which a plasmid/chromosome call from any classifier is not
trusted - the "Ambiguous" tier's second trigger in the four-tier framework
(CLAUDE.md, issue #23): "classifiers disagree, OR contig below the validated
length floor".

Length is read from the same Unicycler `length=` assembly header field
bin/circularity_coverage.py already parses for depth/circularity (see that
script's docstring for the exact header format and why it's the source of
truth here rather than a classifier's own reported contig size column).

Same join-key and sample-scoping approach as bin/classifier_agreement.py
(issue #24, see that file's docstring for why).
"""

import argparse
import gzip
import re

import pandas as pd

__version__ = "0.1.0"

HEADER_RE = re.compile(r"^>(\S+).*?length=(\d+)")


def first_token(contig_id: str) -> str:
    return str(contig_id).split()[0]


def load_hamronization(path: str, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["input_file_name"] == sample_id].copy()
    df["contig_key"] = df["input_sequence_id"].map(first_token)
    return df


def load_contig_lengths(path: str) -> dict:
    """contig_key -> int length, or None if the header has no length= field"""
    opener = gzip.open if path.endswith(".gz") else open
    lengths = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = HEADER_RE.match(line.strip())
            contig_key = first_token(line[1:].strip())
            lengths[contig_key] = int(match.group(2)) if match else None
    return lengths


def main():
    parser = argparse.ArgumentParser(
        prog="contig_length_floor",
        description="Stage 5d: per-ARG validated contig-length floor evidence.",
    )
    parser.add_argument("--hamronization", required=True, help="hAMRonization merged ARG report (tsv, all samples)")
    parser.add_argument("--sample-id", required=True, help="input_file_name value identifying this sample's rows")
    parser.add_argument("--assembly", required=True, help="Assembly fasta (BACASS.out.assembly, may be gzipped)")
    parser.add_argument(
        "--min-length",
        type=int,
        default=1000,
        help=(
            "Contigs shorter than this (bp) are flagged below_length_floor=True. Default 1000bp: a "
            "commonly used floor in short-read plasmid-classification literature below which assembly "
            "fragmentation and misassembly make a plasmid/chromosome call unreliable - not calibrated "
            "against this pipeline's own data yet. Revisit once stage 7 (issue #22) has real validation "
            "data."
        ),
    )
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"contig_length_floor {__version__}")
    args = parser.parse_args()

    hamronization = load_hamronization(args.hamronization, args.sample_id)
    lengths = load_contig_lengths(args.assembly)

    rows = []
    for _, arg_row in hamronization.iterrows():
        contig_key = arg_row["contig_key"]
        length = lengths.get(contig_key)
        below_floor = length is not None and length < args.min_length
        rows.append(
            {
                "gene_symbol": arg_row.get("gene_symbol"),
                "input_sequence_id": arg_row["input_sequence_id"],
                "contig_key": contig_key,
                "contig_length": length,
                "below_length_floor": below_floor,
            }
        )

    columns = ["gene_symbol", "input_sequence_id", "contig_key", "contig_length", "below_length_floor"]
    pd.DataFrame(rows, columns=columns).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
