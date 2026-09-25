#!/usr/bin/env python3

"""
Stage 7e (issue #35): how stable are the tier calls when the input is simulated reads
instead of the closed genome?

DIAGNOSTIC ONLY. The reads are simulated, so this shows how much depth, base quality and
plasmid copy number move the tier calls; it is not, and must never be reported as, an
accuracy figure. Accuracy comes from bin/evaluate_metrics.py on the closed genomes.

Reference: the closed-genome run (one <genome>.tier_resolution.tsv each).
Simulated: the same genome as simulated reads, assembled and run through the whole pipeline
(one <genome>_d<depth>_q<shift>_p<copies>.tier_resolution.tsv per condition, see
assets/validation/simulate_reads.py).

Re-assembly renames contigs, so calls cannot be matched by contig. They are matched by
(genome, gene_symbol): each gene's set of distinct tiers in the simulated run is compared
with its set in the closed-genome run. A gene present in several copies can therefore
match a different number of copies; only the set of tiers is compared.

Per condition and per reference tier (the tier the gene had in the closed-genome run):
    Recovered      gene called at all in the simulated run
    Same_Tier      recovered with exactly the same set of tiers
    Flipped_Class  recovered, but changed between plasmid (High/Moderate) and Chromosomal
                   with no overlap in tiers - the most serious kind of change
    To_Ambiguous   recovered, tiers differ, and it now includes Ambiguous
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd

__version__ = "0.1.0"

SUFFIX = ".tier_resolution.tsv"
# The pipeline appends the assembler to the sample id (e.g. "..._p1-unicycler"); tolerate it.
SAMPLE_RE = re.compile(r"^(?P<genome>.+)_d(?P<depth>\d+)_q(?P<shift>-?\d+)_p(?P<copies>\d+)(?:-[A-Za-z0-9]+)?$")
PLASMID_TIERS = {"High-confidence plasmid", "Moderate-confidence plasmid"}
CHROMOSOMAL = "Chromosomal"
AMBIGUOUS = "Ambiguous"


def tier_sets(paths):
    """({(sample, gene_symbol): frozenset(tiers)}, [sample, ...]) over the given tier files.

    Samples come from the file names, not from the rows: a run that found no ARGs at all
    (the worst outcome of a degraded input) has an empty file and must still count."""
    sets, samples = {}, []
    for path in paths:
        sample = path.name[: -len(SUFFIX)] if path.name.endswith(SUFFIX) else path.stem
        samples.append(sample)
        df = pd.read_csv(path, sep="\t", dtype=str)
        if not {"gene_symbol", "tier"} <= set(df.columns):
            sys.exit(f"Error: {path.name} needs 'gene_symbol' and 'tier' columns")
        for gene, tiers in df.groupby("gene_symbol")["tier"]:
            sets[(sample, gene)] = frozenset(tiers)
    return sets, sorted(samples)


def classify(ref, sim):
    if sim is None:
        return "lost"
    if sim == ref:
        return "same"
    ref_plasmid, ref_chrom = bool(ref & PLASMID_TIERS), CHROMOSOMAL in ref
    sim_plasmid, sim_chrom = bool(sim & PLASMID_TIERS), CHROMOSOMAL in sim
    if (ref_plasmid and not ref_chrom and sim_chrom and not sim_plasmid) or (
        ref_chrom and not ref_plasmid and sim_plasmid and not sim_chrom
    ):
        return "flipped"
    if AMBIGUOUS in sim and AMBIGUOUS not in ref:
        return "to_ambiguous"
    return "changed_other"


def ref_tier_label(ref):
    return next(iter(ref)) if len(ref) == 1 else "Mixed (gene in several tiers)"


def condition_key(sample):
    """Sample id without the assembler suffix, so it can be matched to the samplesheet ID."""
    m = SAMPLE_RE.match(sample)
    return f"{m['genome']}_d{m['depth']}_q{m['shift']}_p{m['copies']}" if m else sample


def build(reference, simulated, sim_samples):
    ref_by_genome = {}
    for (genome, gene), tiers in reference.items():
        ref_by_genome.setdefault(genome, {})[gene] = tiers

    rows = []
    for sample in sim_samples:
        m = SAMPLE_RE.match(sample)
        if not m:
            sys.exit(f"Error: simulated sample '{sample}' is not named <genome>_d<depth>_q<shift>_p<copies>")
        genome = m["genome"]
        if genome not in ref_by_genome:
            sys.exit(f"Error: no closed-genome reference tier file (or it has no ARG calls) for '{genome}'")
        for gene, ref in ref_by_genome[genome].items():
            sim = simulated.get((sample, gene))
            rows.append({
                "genome": genome, "sample": sample, "Depth": int(m["depth"]), "Quality_Shift": int(m["shift"]),
                "Plasmid_Copies": int(m["copies"]), "gene_symbol": gene, "Reference_Tier": ref_tier_label(ref),
                "Reference_Tiers": ";".join(sorted(ref)), "Simulated_Tiers": ";".join(sorted(sim or [])),
                "Outcome": classify(ref, sim),
            })
    return pd.DataFrame(rows)


def summarise(detail):
    keys = ["Depth", "Quality_Shift", "Plasmid_Copies", "Reference_Tier"]
    g = detail.groupby(keys)
    out = g.size().rename("N_Genes").to_frame()
    out["Recovered"] = g["Outcome"].apply(lambda s: int((s != "lost").sum()))
    out["Same_Tier"] = g["Outcome"].apply(lambda s: int((s == "same").sum()))
    out["Flipped_Class"] = g["Outcome"].apply(lambda s: int((s == "flipped").sum()))
    out["To_Ambiguous"] = g["Outcome"].apply(lambda s: int((s == "to_ambiguous").sum()))
    out["Recovered_Pct"] = (100 * out["Recovered"] / out["N_Genes"]).round(1)
    out["Same_Tier_Pct_Of_Recovered"] = (100 * out["Same_Tier"] / out["Recovered"].where(out["Recovered"] > 0)).round(1)
    out["Note"] = "DIAGNOSTIC ONLY - simulated reads, not an accuracy figure"
    return out.reset_index()


def main():
    parser = argparse.ArgumentParser(prog="sensitivity_summary", description="Diagnostic tier-stability summary (stage 7e).")
    parser.add_argument("--reference", required=True, nargs="+", type=Path, help="closed-genome *.tier_resolution.tsv files")
    parser.add_argument("--simulated", required=True, nargs="+", type=Path, help="simulated-read *.tier_resolution.tsv files")
    parser.add_argument("--samplesheet", help="sim_samplesheet.csv; samples in it with no tier file count as every gene lost")
    parser.add_argument("--output", required=True, help="Summary TSV (per condition and reference tier)")
    parser.add_argument("--detail-output", help="Optional per-gene TSV")
    parser.add_argument("--version", action="version", version=f"sensitivity_summary {__version__}")
    args = parser.parse_args()

    reference, _ = tier_sets(args.reference)
    simulated, sim_samples = tier_sets(args.simulated)
    if args.samplesheet:
        # A sample whose assembly failed writes no tier file at all. It must still count
        # (every gene lost), not vanish from the denominator.
        present = {condition_key(s) for s in sim_samples}
        expected = pd.read_csv(args.samplesheet, dtype=str)["ID"]
        missing = [i for i in expected if i not in present]
        sim_samples = sorted(sim_samples + missing)
        if missing:
            print(f"{len(missing)} sample(s) produced no output and count as every gene lost: {missing}")
    detail = build(reference, simulated, sim_samples)
    if detail.empty:
        sys.exit("Error: nothing to compare (no genes in the reference tier files)")
    summarise(detail).to_csv(args.output, sep="\t", index=False)
    if args.detail_output:
        detail.to_csv(args.detail_output, sep="\t", index=False)
    print(f"Compared {detail['sample'].nunique()} simulated sample(s) against "
          f"{detail['genome'].nunique()} closed genome(s) -> {args.output}")


if __name__ == "__main__":
    main()
