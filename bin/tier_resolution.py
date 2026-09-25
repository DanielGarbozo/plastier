#!/usr/bin/env python3

"""
Stage 5f (issue #29): tier-resolution rule engine.

Combines the five signals computed by stages 5a-5e into exactly one of the
four confidence tiers per ARG, per the table in CLAUDE.md / issue #23:

| Tier                        | Rule                                                                              |
|------------------------------|------------------------------------------------------------------------------------|
| High-confidence plasmid      | All 3 classifiers agree + replicon detected + circularised or coverage matches   |
| Moderate-confidence plasmid  | Majority agree + >=1 mobility/replication marker, circularisation/coverage missing |
| Ambiguous                    | Classifiers disagree, or contig below the validated length floor                  |
| Chromosomal                  | Classifiers agree on chromosomal origin, or gene inside a typed SCCmec cassette   |

The table as written does not give an explicit rule for every combination
the five upstream signals can actually produce. Those gaps were filled with
an explicit, documented judgement call rather than left as undefined
behaviour - each is called out below and in docs/decisions.md, since this is
exactly the kind of decision CLAUDE.md requires a decisions.md entry for
(this table is what stage 7's benchmark validates).

Precedence, evaluated top to bottom, first match wins:

1. SCCmec override (stage 5e/#28) -> Chromosomal, unconditionally. If the
   gene physically sits inside a typed SCCmec cassette, that is direct
   structural evidence explaining away any classifier confusion - it wins
   over every other signal, including a short-contig or disagreement result
   that would otherwise route to Ambiguous. [Judgement call: CLAUDE.md's
   table doesn't state precedence between the two "OR" conditions across
   different tiers explicitly; this reads as the more informative signal.]

2. Below the contig-length floor (stage 5d/#27) -> Ambiguous. Evaluated
   before classifier agreement itself, since a too-short contig makes the
   classifier calls it would otherwise be judged on untrustworthy in the
   first place.

3. Classifier disagreement (stage 5a/#24) -> Ambiguous, per the table.

4. Unanimous chromosome -> Chromosomal, per the table.

5. Majority chromosome -> Ambiguous. [Judgement call: the table only defines
   a "moderate" tier for majority-*plasmid* agreement; it has no equivalent
   partial-confidence chromosomal tier. Routing majority-chromosome to
   Chromosomal would claim more certainty than 2/3 classifier agreement
   supports for a category (chromosomal/clonal attribution) the table
   otherwise treats as requiring full agreement - Ambiguous is the
   conservative choice.]

6. Unanimous or majority plasmid, with a mobility/replication marker
   (stage 5b/#25):
     - AND unanimous AND (circularised OR coverage matches) (stage 5c/#26)
       -> High-confidence plasmid (all three of the table's High conditions).
     - otherwise -> Moderate-confidence plasmid (majority + marker, missing
       the circularisation/coverage evidence the table requires for High -
       this also covers a *unanimous* plasmid call that has a marker but
       lacks circularisation/coverage evidence, since unanimous agreement
       trivially satisfies "majority agree").

7. Unanimous or majority plasmid, no marker at all -> Ambiguous. [Judgement
   call: the table requires >=1 marker for *both* plasmid tiers: a plasmid
   classifier signal with zero replicon/relaxase evidence doesn't meet the
   floor for either, so it can't be called Moderate, and isn't classifier
   disagreement either, so it isn't literally what the table calls
   Ambiguous - routed there anyway as the closest fit: real plasmid signal
   that isn't strong enough to act on.]

8. Anything else (e.g. no classifier data at all for this ARG's contig) ->
   Ambiguous, as a safe default.
"""

import argparse

import pandas as pd

__version__ = "0.1.0"

JOIN_KEYS = ["gene_symbol", "input_sequence_id"]


def load(path: str, columns: list) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)
    return df[JOIN_KEYS + columns]


def to_bool(series: pd.Series) -> pd.Series:
    return series.map(lambda v: str(v).strip().lower() == "true")


def resolve_tier(row) -> str:
    if row["sccmec_override"]:
        return "Chromosomal"
    if row["below_length_floor"]:
        return "Ambiguous"

    agreement = row["agreement"]
    if agreement == "disagreement" or agreement not in (
        "unanimous_chromosome",
        "unanimous_plasmid",
        "majority_chromosome",
        "majority_plasmid",
    ):
        return "Ambiguous"
    if agreement == "unanimous_chromosome":
        return "Chromosomal"
    if agreement == "majority_chromosome":
        return "Ambiguous"

    # agreement is unanimous_plasmid or majority_plasmid from here on
    if not row["has_marker"]:
        return "Ambiguous"
    if agreement == "unanimous_plasmid" and (row["circular"] or row["coverage_matches_plasmid"]):
        return "High-confidence plasmid"
    return "Moderate-confidence plasmid"


def main():
    parser = argparse.ArgumentParser(
        prog="tier_resolution",
        description="Stage 5f: resolve the five upstream evidence signals into a four-tier confidence call per ARG.",
    )
    parser.add_argument("--classifier-agreement", required=True)
    parser.add_argument("--mobility-markers", required=True)
    parser.add_argument("--circularity-coverage", required=True)
    parser.add_argument("--contig-length-floor", required=True)
    parser.add_argument("--sccmec-override", required=True)
    parser.add_argument("--output", required=True, help="Output tsv path")
    parser.add_argument("--version", action="version", version=f"tier_resolution {__version__}")
    args = parser.parse_args()

    agreement_df = load(
        args.classifier_agreement,
        ["mobsuite_call", "platon_call", "rfplasmid_call", "n_classifiers", "n_plasmid_votes", "n_chromosome_votes", "agreement"],
    )
    markers_df = load(args.mobility_markers, ["rep_type", "relaxase_type", "predicted_mobility", "has_marker"])
    coverage_df = load(args.circularity_coverage, ["circular", "coverage_ratio", "coverage_matches_plasmid", "high_confidence_evidence"])
    length_df = load(args.contig_length_floor, ["contig_length", "below_length_floor"])
    sccmec_df = load(args.sccmec_override, ["sccmec_override", "sccmec_element_span"])

    merged = agreement_df
    for other in (markers_df, coverage_df, length_df, sccmec_df):
        merged = merged.merge(other, on=JOIN_KEYS, how="left")

    for col in ("has_marker", "circular", "coverage_matches_plasmid", "below_length_floor", "sccmec_override"):
        merged[col] = to_bool(merged[col])

    merged["tier"] = merged.apply(resolve_tier, axis=1)

    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
