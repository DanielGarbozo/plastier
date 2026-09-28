#!/usr/bin/env python3

"""
Download the closed reference genomes in closed_genome_ground_truth.csv and write
the samplesheet that `--assemblies` takes (columns: sample,fasta).

Each genome is fetched from NCBI Datasets as its own RefSeq FASTA and saved as
<outdir>/<Genome_Accession>.fna.gz. Every plasmid accession the ground truth lists
for that genome must appear as a sequence id in the FASTA - stage 7 scores ARGs by
exactly that id (bin/evaluate_metrics.py), so a mismatch here would silently
corrupt the benchmark and is treated as an error.

Headers are rewritten to Unicycler's style, ">ACCESSION length=N depth=1.00x circular=true",
(sequences untouched, accession still the first token). Stage 5's circularity evidence
(bin/circularity_coverage.py) reads exactly those fields, and an NCBI header has none of
them, which would make the High-confidence plasmid tier unreachable in the benchmark.
These are closed, complete genomes, so circular=true is a fact - and it is set on EVERY
replicon, chromosome included, so it says nothing about which are plasmids. depth=1.00x
everywhere means no coverage evidence is provided (a real plasmid's raised depth is what
the simulated-read analysis, #35, varies). Pass --keep-ncbi-headers to skip the rewrite.

    python3 assets/validation/fetch_reference_genomes.py \\
        --outdir /scratch/$USER/plastier_refs --samplesheet refs.csv

Already-downloaded genomes are reused, so re-running resumes an interrupted fetch.
"""

import argparse
import csv
import gzip
import io
import re
import sys
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

DATASETS_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
    "?include_annotation_type=GENOME_FASTA&filename={acc}.zip"
)
REQUEST_DELAY_S = 0.5  # stay well under NCBI's unauthenticated 5 requests/s
RETRIES = 3


def read_ground_truth(path):
    plasmids = {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh):
            plasmids.setdefault(row["Genome_Accession"], set()).add(row["Plasmid_Accession"])
    return plasmids


def download_fasta(acc):
    url = DATASETS_URL.format(acc=acc)
    for attempt in range(1, RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=120) as resp:
                archive = zipfile.ZipFile(io.BytesIO(resp.read()))
            names = [n for n in archive.namelist() if n.endswith("_genomic.fna")]
            if len(names) != 1:
                sys.exit(f"Error: expected one *_genomic.fna in the {acc} download, found {names}")
            return archive.read(names[0])
        except (urllib.error.URLError, zipfile.BadZipFile, TimeoutError) as err:
            if attempt == RETRIES:
                sys.exit(f"Error: could not download {acc} after {RETRIES} attempts: {err}")
            time.sleep(2 * attempt)



def unicycler_headers(fasta_bytes):
    """Rewrite each header to '>ACC length=N depth=1.00x circular=true'; idempotent."""
    records, name, seq = [], None, []
    for line in fasta_bytes.splitlines():
        if line.startswith(b">"):
            if name is not None:
                records.append((name, seq))
            name, seq = line[1:].split()[0], []
        elif name is not None:
            seq.append(line)
    if name is not None:
        records.append((name, seq))
    out = []
    for name, seq in records:
        length = sum(len(part) for part in seq)
        out.append(b">" + name + b" length=" + str(length).encode() + b" depth=1.00x circular=true")
        out.extend(seq)
    return b"\n".join(out) + b"\n"


def sequence_ids(fasta_bytes):
    return {line[1:].split()[0].decode() for line in fasta_bytes.splitlines() if line.startswith(b">")}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ground-truth", default=str(Path(__file__).with_name("closed_genome_ground_truth.csv")))
    parser.add_argument("--outdir", required=True, type=Path, help="Directory for the <accession>.fna.gz files")
    parser.add_argument("--samplesheet", required=True, type=Path, help="Output csv for --assemblies")
    parser.add_argument("--only", nargs="+", help="Restrict to these genome accessions (for testing)")
    parser.add_argument("--keep-ncbi-headers", action="store_true",
                        help="Do not rewrite headers to Unicycler style (High-confidence tier then cannot occur)")
    args = parser.parse_args()

    plasmids = read_ground_truth(args.ground_truth)
    genomes = sorted(args.only) if args.only else sorted(plasmids)
    unknown = [g for g in genomes if g not in plasmids]
    if unknown:
        sys.exit(f"Error: not in the ground truth: {unknown}")

    args.outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for i, acc in enumerate(genomes, 1):
        target = args.outdir / f"{acc}.fna.gz"
        if target.exists():
            with gzip.open(target, "rb") as fh:
                fasta = fh.read()
            status = "cached"
        else:
            time.sleep(REQUEST_DELAY_S)
            fasta = download_fasta(acc)
            status = "downloaded"

        ids = sequence_ids(fasta)
        missing = plasmids[acc] - ids
        if missing:
            target.unlink(missing_ok=True)
            sys.exit(f"Error: {acc}: plasmid accession(s) {sorted(missing)} not among the FASTA's sequence ids {sorted(ids)}")
        annotated = fasta if args.keep_ncbi_headers else unicycler_headers(fasta)
        if status == "downloaded" or annotated != fasta:
            with gzip.open(target, "wb") as fh:
                fh.write(annotated)
            if status == "cached":
                status = "cached, headers rewritten"

        print(f"[{i}/{len(genomes)}] {acc}: {status}, {len(ids)} replicons ({len(plasmids[acc])} plasmid)")
        rows.append((acc, str(target.resolve())))

    with open(args.samplesheet, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "fasta"])
        writer.writerows(rows)
    print(f"Wrote {len(rows)} genomes to {args.samplesheet}")


if __name__ == "__main__":
    main()
