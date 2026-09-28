#!/usr/bin/env python3

"""
Download the Phase 3 candidate runs from ENA and write a samplesheet for `--input`.

Why this exists instead of `--sra_ids` (stage 1, nf-core/fetchngs): on the Slurm cluster this
was developed on, DNS does not resolve inside Apptainer containers (the images have no
/etc/resolv.conf to bind onto), so fetchngs' wget-in-a-container download step fails for every
sample. Downloading on the host avoids that, and every file is checked against the MD5 that ENA
publishes, which the `--sra_ids` path does too. Stage 2 onwards is unaffected.

    python3 assets/phase3/fetch_reads.py --ids assets/phase3/sra_ids_candidates.txt \\
        --outdir /scratch/$USER/plastier_phase3_reads --samplesheet <outdir>/samplesheet.csv

Already-downloaded, verified files are reused, so an interrupted fetch resumes.
"""

import argparse
import csv
import hashlib
import io
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

ENA = "https://www.ebi.ac.uk/ena/portal/api/filereport"
RETRIES = 3


def filereport(run):
    query = urllib.parse.urlencode({
        "accession": run, "result": "read_run", "format": "tsv",
        "fields": "run_accession,scientific_name,library_layout,fastq_ftp,fastq_md5",
    })
    with urllib.request.urlopen(f"{ENA}?{query}", timeout=60) as resp:
        rows = list(csv.DictReader(io.StringIO(resp.read().decode()), delimiter="\t"))
    if len(rows) != 1:
        sys.exit(f"Error: ENA returned {len(rows)} records for {run}")
    return rows[0]


def md5sum(path):
    digest = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url, target, expected_md5):
    if target.exists() and md5sum(target) == expected_md5:
        return "cached"
    for attempt in range(1, RETRIES + 1):
        try:
            urllib.request.urlretrieve(url, target)
            if md5sum(target) == expected_md5:
                return "downloaded"
            target.unlink()
            raise OSError("md5 mismatch")
        except OSError as err:
            if attempt == RETRIES:
                sys.exit(f"Error: {url}: {err} after {RETRIES} attempts")
            time.sleep(2 * attempt)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ids", required=True, type=Path, help="one run accession per line")
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--samplesheet", required=True, type=Path)
    args = parser.parse_args()

    runs = [line.strip() for line in args.ids.read_text().splitlines() if line.strip()]
    args.outdir.mkdir(parents=True, exist_ok=True)
    sheet = []
    for i, run in enumerate(runs, 1):
        info = filereport(run)
        if info["scientific_name"] != "Staphylococcus aureus" or info["library_layout"] != "PAIRED":
            sys.exit(f"Error: {run} is '{info['scientific_name']}', {info['library_layout']} - not paired S. aureus")
        urls, md5s = info["fastq_ftp"].split(";"), info["fastq_md5"].split(";")
        if len(urls) != 2:
            sys.exit(f"Error: {run} has {len(urls)} FASTQ files, expected 2 (R1/R2 only)")
        paths, status = [], []
        for url, md5 in sorted(zip(urls, md5s)):  # _1 before _2
            target = args.outdir / Path(url).name
            status.append(download("https://" + url, target, md5))
            paths.append(target.resolve())
        print(f"[{i}/{len(runs)}] {run}: {'/'.join(status)}")
        sheet.append((run, paths[0], paths[1]))

    with open(args.samplesheet, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"])
        writer.writerows((run, r1, r2, "NA", "NA", "2.8m") for run, r1, r2 in sheet)
    print(f"Wrote {len(sheet)} samples to {args.samplesheet}")


if __name__ == "__main__":
    main()
