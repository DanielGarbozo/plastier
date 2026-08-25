#!/usr/bin/env python3

"""
Stage 5a (issue #24): classifier-agreement signal.

Intersects stage 3's ARG calls (hAMRonization's merged report) with stage 4's
per-contig plasmid/chromosome calls from MOB-suite, Platon, and RFPlasmid, and
reports, per ARG, how many of the three classifiers called its contig
"plasmid" vs "chromosome". This is the foundational signal stage 5's tier
table (see CLAUDE.md and issue #23) reads first: "all 3 agree", "majority
agree", or "classifiers disagree".

Join key: each tool's contig identifier, truncated to its first
whitespace-delimited token. hAMRonization's ARG callers (ABRicate,
AMRFinderPlus, ...) report only that first token as `input_sequence_id`, but
MOB-suite and RFPlasmid preserve the full FASTA header (id + description) -
joining on the untruncated string would silently match nothing. All three
tools run against the same BACASS assembly, so the first token is a stable,
unambiguous join key across all of them.

Platon's TSV lists only the contigs it calls plasmid (confirmed from
upstream's README, not yet exercised against a real run in this repo - its
~2.8 GB database makes that expensive, see subworkflows/local/
plasmid_classification/main.nf) - so a contig absent from that file is an
implicit chromosome call. Platon input is optional: pass --platon-tsv only
when the run didn't skip it (--plasmid_skip_platon).

Sample scoping: ARG's hamronization report (subworkflows/local/funcscan_arg)
is merged across every sample in the run in one shot (HAMRONIZATION_SUMMARIZE
collects all samples before writing it) - unlike MOB-suite/Platon/RFPlasmid,
whose outputs stay one-file-per-sample all the way through
subworkflows/local/plasmid_classification. Rather than reshape that already-
shipped stage 3 behaviour, this script takes the whole merged report plus
--sample-id and filters to that sample's rows (matched against
`input_file_name`) itself. The calling subworkflow broadcasts the same report
path to every sample's invocation.
"""

import argparse
import sys

import pandas as pd

__version__ = "0.1.0"


def first_token(contig_id: str) -> str:
    return str(contig_id).split()[0]


def load_hamronization(path: str, sample_id: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    df = df[df["input_file_name"] == sample_id].copy()
    df["contig_key"] = df["input_sequence_id"].map(first_token)
    return df


def load_mobsuite(path: str) -> dict:
    """contig_key -> 'plasmid' | 'chromosome'"""
    df = pd.read_csv(path, sep="\t", dtype=str)
    df["contig_key"] = df["contig_id"].map(first_token)
    # molecule_type is one of 'chromosome', 'plasmid', or 'unknown' upstream -
    # anything not literally 'chromosome' is treated as a plasmid call.
    calls = df.set_index("contig_key")["molecule_type"].map(
        lambda v: "chromosome" if str(v).strip().lower() == "chromosome" else "plasmid"
    )
    return calls.to_dict()


def load_rfplasmid(path: str) -> dict:
    """contig_key -> 'plasmid' | 'chromosome'"""
    df = pd.read_csv(path, dtype=str)
    df["contig_key"] = df["contigID"].map(first_token)
    calls = df.set_index("contig_key")["prediction"].map(
        lambda v: "plasmid" if str(v).strip().lower().startswith("p") else "chromosome"
    )
    return calls.to_dict()


def load_platon(path: str | None, chromosome_contig_keys: set) -> dict:
    """contig_key -> 'plasmid' | 'chromosome', or {} if Platon was skipped."""
    if not path:
        return {}
    df = pd.read_csv(path, sep="\t", dtype=str)
    if df.empty:
        # Platon found no plasmid contigs at all - every contig it saw is
        # chromosome, but we only know which contigs it *saw* from the other
        # tools' output, since an empty Platon tsv carries no contig IDs.
        return {key: "chromosome" for key in chromosome_contig_keys}
    id_col = df.columns[0]
    plasmid_keys = set(df[id_col].map(first_token))
    return {key: ("plasmid" if key in plasmid_keys else "chromosome") for key in chromosome_contig_keys}


def agreement_level(calls: list) -> str:
    """calls: list of 'plasmid'/'chromosome' strings from whichever classifiers ran."""
    if not calls:
        return "no_classifier_data"
    n_plasmid = calls.count("plasmid")
    n_chromosome = calls.count("chromosome")
    if n_plasmid == len(calls):
        return "unanimous_plasmid"
    if n_chromosome == len(calls):
        return "unanimous_chromosome"
    if n_plasmid > n_chromosome:
        return "majority_plasmid"
    if n_chromosome > n_plasmid:
        return "majority_chromosome"
    return "disagreement"


def main():
    parser = argparse.ArgumentParser(
        prog="classifier_agreement",
        description="Stage 5a: per-ARG plasmid/chromosome classifier-agreement signal.",
    )
    parser.add_argument("--hamronization", required=True, help="hAMRonization merged ARG report (tsv, all samples)")
    parser.add_argument("--sample-id", required=True, help="input_file_name value identifying this sample's rows")
    parser.add_argument("--mobsuite", required=True, help="MOB-suite contig_report.txt")
    parser.add_argument("--rfplasmid", required=True, help="RFPlasmid prediction.csv")
    parser.add_argument("--platon", default=None, help="Platon tsv (omit if --plasmid_skip_platon was set)")
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"classifier_agreement {__version__}")
    args = parser.parse_args()

    hamronization = load_hamronization(args.hamronization, args.sample_id)
    mobsuite_calls = load_mobsuite(args.mobsuite)
    rfplasmid_calls = load_rfplasmid(args.rfplasmid)
    # Platon only lists plasmid hits, so its universe of "seen" contigs has to
    # be inferred from a tool that reports every contig.
    platon_calls = load_platon(args.platon, set(mobsuite_calls.keys()) | set(rfplasmid_calls.keys()))

    rows = []
    unmatched_contigs = set()
    for _, arg_row in hamronization.iterrows():
        contig_key = arg_row["contig_key"]
        calls = []
        mobsuite_call = mobsuite_calls.get(contig_key)
        rfplasmid_call = rfplasmid_calls.get(contig_key)
        platon_call = platon_calls.get(contig_key)
        for call in (mobsuite_call, platon_call, rfplasmid_call):
            if call is not None:
                calls.append(call)
        if not calls:
            unmatched_contigs.add(contig_key)

        rows.append(
            {
                "gene_symbol": arg_row.get("gene_symbol"),
                "input_sequence_id": arg_row["input_sequence_id"],
                "contig_key": contig_key,
                "mobsuite_call": mobsuite_call,
                "platon_call": platon_call,
                "rfplasmid_call": rfplasmid_call,
                "n_classifiers": len(calls),
                "n_plasmid_votes": calls.count("plasmid"),
                "n_chromosome_votes": calls.count("chromosome"),
                "agreement": agreement_level(calls),
            }
        )

    if unmatched_contigs:
        print(
            f"WARNING: {len(unmatched_contigs)} contig(s) referenced by an ARG call had no match in any "
            f"classifier output (checked by first whitespace-delimited token): {sorted(unmatched_contigs)}",
            file=sys.stderr,
        )

    columns = [
        "gene_symbol",
        "input_sequence_id",
        "contig_key",
        "mobsuite_call",
        "platon_call",
        "rfplasmid_call",
        "n_classifiers",
        "n_plasmid_votes",
        "n_chromosome_votes",
        "agreement",
    ]
    pd.DataFrame(rows, columns=columns).to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
