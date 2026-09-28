#!/usr/bin/env python3

"""
Phase 3 candidate selection: ~20 public S. aureus Illumina runs for the scale-up run.

Queries ENA's portal API for S. aureus paired-end Illumina WGS runs and picks a
diverse, well-behaved subset with a fixed seed. ENA keeps growing, so re-running this
later can return different runs; the committed sra_ids_candidates.txt and
candidates_metadata.tsv are the record of what was actually selected and reviewed.

Selection rules (all applied here, none by hand):
  1. scientific_name is exactly "Staphylococcus aureus", library_strategy WGS,
     library_source GENOMIC, platform ILLUMINA, layout PAIRED.
  2. 250-450 Mb of sequence (~85-150x of a ~2.9 Mb genome) and 120-151 bp mean read
     length: enough depth to assemble plasmids, small enough for a 20-genome run.
  3. country, collection_date, host and isolation_source are all present, so every
     candidate can be described and discussed.
  4. Not one of the Phase 0 pilot runs, and not from a study a pilot run came from.
  5. --human human-host runs plus one run for each of the animal hosts in ANIMAL_HOSTS.
     At most one run per study within the human runs, and the human runs avoid every
     study an animal run came from. Some animal hosts have only one or two studies in
     ENA (several hosts can share one), so an animal run prefers a study not used yet
     but may fall back to any study of that host. Countries are spread first (none
     repeats until every available country has been used); studies are visited in
     seeded random order.

    python3 assets/phase3/select_candidates.py --outdir assets/phase3
"""

import argparse
import io
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ENA = "https://www.ebi.ac.uk/ena/portal/api/search"
QUERY = (
    'tax_eq(1280) AND library_strategy="WGS" AND library_source="GENOMIC" '
    'AND instrument_platform="ILLUMINA" AND library_layout="PAIRED"'
)
FIELDS = (
    "run_accession,scientific_name,study_accession,sample_accession,center_name,instrument_model,"
    "read_count,base_count,country,collection_date,first_public,fastq_bytes,isolation_source,host"
)
ANIMAL_HOSTS = ["Felis catus", "Bos taurus", "Canis lupus familiaris", "Equus caballus", "Sus scrofa domesticus"]
SEED = 42


def fetch(params):
    url = ENA + "?" + urllib.parse.urlencode(params)
    with urllib.request.urlopen(url, timeout=300) as resp:
        return pd.read_csv(io.BytesIO(resp.read()), sep="\t", dtype=str)


def pick_spread(pool, n, rng_seed):
    """One run per study; use each country once before any country twice."""
    pool = pool.sample(frac=1, random_state=rng_seed).drop_duplicates("study_accession")
    chosen, used_countries = [], set()
    while len(chosen) < n and len(pool):
        fresh = pool[~pool["country_main"].isin(used_countries)]
        row = (fresh if len(fresh) else pool).iloc[0]
        chosen.append(row)
        used_countries.add(row["country_main"])
        pool = pool[pool["study_accession"] != row["study_accession"]]
    return pd.DataFrame(chosen)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--outdir", type=Path, default=Path(__file__).parent)
    parser.add_argument("--pilot", type=Path, default=Path(__file__).parents[2] / "pilot_sra.csv",
                        help="Phase 0 pilot run accessions, one per line (excluded, with their studies)")
    parser.add_argument("--human", type=int, default=15, help="number of human-host runs")
    args = parser.parse_args()

    runs = fetch({"result": "read_run", "query": QUERY, "fields": FIELDS, "format": "tsv", "limit": 0})
    print(f"ENA returned {len(runs)} S. aureus paired Illumina WGS runs")

    pilot_ids = [line.strip() for line in args.pilot.read_text().splitlines() if line.strip()]
    pilot_studies = set(runs.loc[runs["run_accession"].isin(pilot_ids), "study_accession"])

    runs["bases"] = pd.to_numeric(runs["base_count"], errors="coerce")
    runs["read_len"] = runs["bases"] / (pd.to_numeric(runs["read_count"], errors="coerce") * 2)
    runs["country_main"] = runs["country"].str.split(":").str[0].str.strip()
    ok = runs[
        (runs["scientific_name"] == "Staphylococcus aureus")
        & runs["bases"].between(2.5e8, 4.5e8)
        & runs["read_len"].between(120, 151)
        & runs[["country_main", "collection_date", "host", "isolation_source"]].notna().all(axis=1)
        & ~runs["run_accession"].isin(pilot_ids)
        & ~runs["study_accession"].isin(pilot_studies)
    ]
    print(f"{len(ok)} runs pass the filters ({ok['study_accession'].nunique()} studies)")

    # Animal hosts first (fewest studies), then humans from studies not yet used.
    animals, used_studies = [], set()
    for host in ANIMAL_HOSTS:
        host_runs = ok[ok["host"] == host]
        fresh = host_runs[~host_runs["study_accession"].isin(used_studies)]
        pick = pick_spread(fresh if len(fresh) else host_runs, 1, SEED)
        animals.append(pick)
        used_studies |= set(pick["study_accession"])
    human = pick_spread(ok[(ok["host"] == "Homo sapiens") & ~ok["study_accession"].isin(used_studies)], args.human, SEED)
    chosen = pd.concat([human, *animals], ignore_index=True)
    assert chosen["run_accession"].is_unique and human["study_accession"].is_unique
    assert not set(human["study_accession"]) & used_studies

    cols = ["run_accession", "study_accession", "sample_accession", "host", "country", "collection_date",
            "isolation_source", "center_name", "instrument_model", "read_count", "base_count", "fastq_bytes"]
    args.outdir.mkdir(parents=True, exist_ok=True)
    chosen[cols].to_csv(args.outdir / "candidates_metadata.tsv", sep="\t", index=False)
    (args.outdir / "sra_ids_candidates.txt").write_text("\n".join(chosen["run_accession"]) + "\n")
    print(chosen[["run_accession", "host", "country", "collection_date", "isolation_source"]].to_string(index=False))


if __name__ == "__main__":
    main()
