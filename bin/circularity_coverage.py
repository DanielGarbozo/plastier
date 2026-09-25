#!/usr/bin/env python3

"""
Stage 5c (issue #26): circularisation / coverage-ratio evidence.

Reports, per ARG, whether its contig is circularised and whether its
coverage ratio matches a multi-copy plasmid rather than baseline chromosomal
depth. This is what distinguishes "High-confidence plasmid" from
"Moderate-confidence plasmid" in the four-tier framework (CLAUDE.md, issue
#23) - both tiers require classifier agreement and a mobility/replication
marker, but only the high-confidence tier also requires circularisation or a
matching coverage ratio.

Sourced directly from the assembly FASTA's own contig headers, not from
MOB-suite's contig_report.txt (whose circularity_status is frequently "not
tested" for short-read-only assemblies) or Platon (frequently skipped - see
subworkflows/local/plasmid_classification/main.nf - and its coverage column
was never exercised against a real run in this repo). Unicycler (this
pipeline's default short-read assembler) writes headers like
">1 length=384426 depth=0.99x" or ">3 length=5320 depth=3.12x circular=true"
- confirmed against a real assembly produced by this repo's own bacass test
fixtures (subworkflows/local/bacass), not just documentation. Unicycler
already normalises depth to the assembly's median contig depth, so the
chromosome sits near 1.00x and a multi-copy plasmid reads measurably higher -
the raw depth value from the header *is* the coverage ratio, no extra
computation needed.

Only Unicycler's header convention is handled - an assembly produced by any
other bacass assembler path has no depth=/circular= annotations, so every
contig comes back with coverage_ratio=None, circular=False (missing
evidence, not a false negative): see --assembler in subworkflows/local/bacass.

Same join-key and sample-scoping approach as bin/classifier_agreement.py
(issue #24, see that file's docstring for why).
"""

import argparse
import gzip
import re

import pandas as pd

__version__ = "0.1.0"

HEADER_RE = re.compile(r"^>(\S+).*?depth=([\d.]+)x(.*)$")


def first_token(contig_id: str) -> str:
    return str(contig_id).split()[0]


def load_hamronization(path: str, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["input_file_name"] == sample_id].copy()
    df["contig_key"] = df["input_sequence_id"].map(first_token)
    return df


def load_assembly_headers(path: str) -> dict:
    """contig_key -> {coverage_ratio: float|None, circular: bool}"""
    opener = gzip.open if path.endswith(".gz") else open
    headers = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            match = HEADER_RE.match(line.strip())
            if not match:
                # Not a Unicycler-style header (different assembler) - no
                # depth/circular evidence available for this contig.
                contig_key = first_token(line[1:].strip())
                headers[contig_key] = {"coverage_ratio": None, "circular": False}
                continue
            contig_id, depth, rest = match.groups()
            headers[first_token(contig_id)] = {
                "coverage_ratio": float(depth),
                "circular": "circular=true" in rest,
            }
    return headers


def main():
    parser = argparse.ArgumentParser(
        prog="circularity_coverage",
        description="Stage 5c: per-ARG circularisation / coverage-ratio evidence from the assembly headers.",
    )
    parser.add_argument("--hamronization", required=True, help="hAMRonization merged ARG report (tsv, all samples)")
    parser.add_argument("--sample-id", required=True, help="input_file_name value identifying this sample's rows")
    parser.add_argument("--assembly", required=True, help="Assembly fasta (BACASS.out.assembly, may be gzipped)")
    parser.add_argument(
        "--coverage-ratio-threshold",
        type=float,
        default=1.5,
        help=(
            "Minimum depth= value (relative to the assembly's own median, per Unicycler's normalisation) "
            "to count as plasmid-like multi-copy coverage rather than baseline chromosomal depth. "
            "Default 1.5x: real chromosomal contigs vary roughly 0.9-1.1x in practice (see the fixture "
            "evidence cited in this script's docstring), so 1.5x is a conservative floor above that noise "
            "band, not a literal copy-number claim. Revisit once stage 7 (issue #22) has real validation "
            "data to calibrate against."
        ),
    )
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"circularity_coverage {__version__}")
    args = parser.parse_args()

    hamronization = load_hamronization(args.hamronization, args.sample_id)
    headers = load_assembly_headers(args.assembly)

    rows = []
    for _, arg_row in hamronization.iterrows():
        contig_key = arg_row["contig_key"]
        info = headers.get(contig_key, {"coverage_ratio": None, "circular": False})
        coverage_ratio = info["coverage_ratio"]
        coverage_matches_plasmid = coverage_ratio is not None and coverage_ratio >= args.coverage_ratio_threshold
        rows.append(
            {
                "gene_symbol": arg_row.get("gene_symbol"),
                "input_sequence_id": arg_row["input_sequence_id"],
                "contig_key": contig_key,
                "circular": info["circular"],
                "coverage_ratio": coverage_ratio,
                "coverage_matches_plasmid": coverage_matches_plasmid,
                "high_confidence_evidence": info["circular"] or coverage_matches_plasmid,
            }
        )

    columns = [
        "gene_symbol",
        "input_sequence_id",
        "contig_key",
        "circular",
        "coverage_ratio",
        "coverage_matches_plasmid",
        "high_confidence_evidence",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
