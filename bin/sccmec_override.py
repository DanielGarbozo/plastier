#!/usr/bin/env python3

"""
Stage 5e (issue #28): SCCmec-cassette override to the chromosomal tier.

Reports, per ARG, whether its coordinates fall inside a typed SCCmec
cassette on the same contig - the four-tier framework's (CLAUDE.md, issue
#23) route to the "Chromosomal" tier regardless of what the plasmid
classifiers said: SCCmec carries its own integrase/mobility genes, so a
naive plasmid classifier can mistake an ARG sitting inside it for mobile
signal when it is actually chromosomal and clonally inherited.

Deliberately uses SCCmecExtractor (subworkflows/local/typing), not
staphopia-sccmec (also in that subworkflow, used for typing elsewhere) -
staphopia-sccmec only reports whether the genome has each SCCmec type
present/absent, confirmed by inspecting its real --json output, which
carries no contig or coordinate information at all. SCCmecExtractor's
sccmec_unified_report.tsv gives the actual att-site boundaries per contig
(Contig, AttL_Start/End, AttR_Start/End columns) - the cassette span is
taken as [min, max] of those four coordinates per row, which is robust to
either attachment-site orientation without needing to interpret which end
is "left" vs "right" in genomic terms.

ARG coordinates come from hAMRonization's input_gene_start/input_gene_stop
columns. Overlap is a standard closed-interval test: the ARG and the
cassette span intersect if the ARG doesn't end before the cassette starts
and doesn't start after the cassette ends.

Same join-key and sample-scoping approach as bin/classifier_agreement.py
(issue #24, see that file's docstring for why).
"""

import argparse

import pandas as pd

__version__ = "0.1.0"


def first_token(contig_id: str) -> str:
    return str(contig_id).split()[0]


def load_hamronization(path: str, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["input_file_name"] == sample_id].copy()
    df["contig_key"] = df["input_sequence_id"].map(first_token)
    return df


def load_sccmec_elements(path: str) -> dict:
    """contig_key -> list of (start, end) cassette spans (extracted elements only)"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["Status"] == "extracted"].copy()
    elements = {}
    coord_cols = ["AttL_Start", "AttL_End", "AttR_Start", "AttR_End"]
    for _, row in df.iterrows():
        contig_key = first_token(row["Contig"])
        coords = [int(row[c]) for c in coord_cols if str(row[c]).strip() not in ("", "-", "nan")]
        if not coords:
            continue
        elements.setdefault(contig_key, []).append((min(coords), max(coords)))
    return elements


def overlaps(arg_start, arg_end, element_start: int, element_end: int) -> bool:
    return arg_start <= element_end and arg_end >= element_start


def main():
    parser = argparse.ArgumentParser(
        prog="sccmec_override",
        description="Stage 5e: per-ARG SCCmec-cassette override to the chromosomal tier.",
    )
    parser.add_argument("--hamronization", required=True, help="hAMRonization merged ARG report (tsv, all samples)")
    parser.add_argument("--sample-id", required=True, help="input_file_name value identifying this sample's rows")
    parser.add_argument("--sccmecextractor", required=True, help="SCCmecExtractor sccmec_unified_report.tsv")
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"sccmec_override {__version__}")
    args = parser.parse_args()

    hamronization = load_hamronization(args.hamronization, args.sample_id)
    elements = load_sccmec_elements(args.sccmecextractor)

    rows = []
    for _, arg_row in hamronization.iterrows():
        contig_key = arg_row["contig_key"]
        sccmec_override = False
        element_span = None

        gene_start_raw = arg_row.get("input_gene_start")
        gene_stop_raw = arg_row.get("input_gene_stop")
        has_coords = (
            pd.notna(gene_start_raw)
            and pd.notna(gene_stop_raw)
            and str(gene_start_raw).strip() != ""
            and str(gene_stop_raw).strip() != ""
        )

        if has_coords and contig_key in elements:
            arg_start, arg_end = sorted((int(float(gene_start_raw)), int(float(gene_stop_raw))))
            for element_start, element_end in elements[contig_key]:
                if overlaps(arg_start, arg_end, element_start, element_end):
                    sccmec_override = True
                    element_span = f"{element_start}-{element_end}"
                    break

        rows.append(
            {
                "gene_symbol": arg_row.get("gene_symbol"),
                "input_sequence_id": arg_row["input_sequence_id"],
                "contig_key": contig_key,
                "sccmec_override": sccmec_override,
                "sccmec_element_span": element_span,
            }
        )

    columns = ["gene_symbol", "input_sequence_id", "contig_key", "sccmec_override", "sccmec_element_span"]
    pd.DataFrame(rows, columns=columns).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
