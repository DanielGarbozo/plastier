#!/usr/bin/env python3

"""
Stage 7e (issue #35): simulate Illumina reads from the closed reference genomes, for the
sensitivity analysis. DIAGNOSTIC ONLY - the reads are simulated, so nothing derived from
them is an accuracy figure (see docs/stage7_validation.md).

For every genome x condition it writes <outdir>/<sample>_R1.fastq.gz / _R2.fastq.gz
(paired-end, 2x150 bp, ART HiSeq 2500 error profile) and one line of the samplesheet
that the pipeline takes with --input.

A condition is depth, quality shift and plasmid copy number, written
    d<depth>_q<shift>_p<copy_number>       e.g. d20_q0_p1   d50_q-15_p3
- depth: mean chromosome coverage in x.
- quality shift: added to every base quality (ART -qs/-qs2); q-15 lowers HiSeq's Q~35 to
  about Q20, i.e. roughly 1% error, to mimic a poor run.
- plasmid copy number: plasmid replicons are simulated at depth x copy_number. Real
  plasmids are usually present in more than one copy per chromosome, and stage 5's
  coverage-ratio evidence reads exactly that, so this is an assumption of the simulation,
  not a measured value: p1 (same depth as the chromosome) is the case where that
  evidence is absent; p3 is a plausible multi-copy plasmid. Compare them, do not pick one.

Replicons are simulated one at a time (ART's -f applies to a whole file), and the plasmids
are identified by the ground-truth accessions, so the split needs no guesswork.

    python3 assets/validation/simulate_reads.py \\
        --assemblies <refs>/assemblies.csv --outdir <simdir> --genomes GCF_... GCF_... \\
        --conditions d20_q0_p1 d50_q0_p1 d20_q-15_p1 d50_q0_p3
"""

import argparse
import csv
import gzip
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

ART_IMAGE_URL = "https://depot.galaxyproject.org/singularity/art:2016.06.05--h589041f_9"
CONDITION_RE = re.compile(r"^d(\d+)_q(-?\d+)_p(\d+)$")
READ_LEN, FRAG_MEAN, FRAG_SD = 150, 350, 30
SEED = 1234


def parse_condition(text):
    m = CONDITION_RE.match(text)
    if not m:
        sys.exit(f"Error: condition '{text}' is not of the form d<depth>_q<shift>_p<copies>, e.g. d20_q0_p1")
    return int(m[1]), int(m[2]), int(m[3])


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    name, seq = None, []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name:
                    yield name, "".join(seq)
                name, seq = line[1:].split()[0], []
            else:
                seq.append(line.strip())
    if name:
        yield name, "".join(seq)


def plasmid_accessions(ground_truth):
    plasmids = {}
    with open(ground_truth, newline="") as fh:
        for row in csv.DictReader(fh):
            plasmids.setdefault(row["Genome_Accession"], set()).add(row["Plasmid_Accession"])
    return plasmids


def art(image, workdir, fasta, prefix, depth, shift):
    cmd = ["apptainer", "exec", "--bind", str(workdir), str(image), "art_illumina",
           "-ss", "HS25", "-i", str(fasta), "-p", "-l", str(READ_LEN), "-m", str(FRAG_MEAN), "-s", str(FRAG_SD),
           "-f", str(depth), "-qs", str(shift), "-qs2", str(shift), "-rs", str(SEED), "-na", "-o", str(prefix)]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode:
        sys.exit(f"Error: ART failed:\n{proc.stderr[-800:]}")


def simulate(image, workdir, outdir, genome, fasta, plasmids, condition):
    depth, shift, copies = parse_condition(condition)
    sample = f"{genome}_{condition}"
    r1, r2 = outdir / f"{sample}_R1.fastq.gz", outdir / f"{sample}_R2.fastq.gz"
    if r1.exists() and r2.exists():
        return sample, r1, r2
    with gzip.open(r1, "wt") as out1, gzip.open(r2, "wt") as out2:
        for name, seq in read_fasta(fasta):
            rep_depth = depth * (copies if name in plasmids else 1)
            rep_fa = workdir / f"{sample}.{name}.fa"
            rep_fa.write_text(f">{name}\n{seq}\n")
            prefix = workdir / f"{sample}.{name}."
            art(image, workdir, rep_fa, prefix, rep_depth, shift)
            for src, out in ((f"{prefix}1.fq", out1), (f"{prefix}2.fq", out2)):
                out.write(Path(src).read_text())
                Path(src).unlink()
            rep_fa.unlink()
            Path(f"{prefix}1.aln").unlink(missing_ok=True)
            Path(f"{prefix}2.aln").unlink(missing_ok=True)
    return sample, r1, r2


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--assemblies", required=True, type=Path, help="assemblies.csv from fetch_reference_genomes.py")
    parser.add_argument("--ground-truth", type=Path, default=Path(__file__).with_name("closed_genome_ground_truth.csv"))
    parser.add_argument("--outdir", required=True, type=Path)
    parser.add_argument("--genomes", required=True, nargs="+", help="Genome accessions to simulate")
    parser.add_argument("--conditions", required=True, nargs="+", help="d<depth>_q<shift>_p<copies> ...")
    parser.add_argument("--art-image", type=Path, help="ART Singularity image (downloaded to <outdir> if omitted)")
    args = parser.parse_args()

    for c in args.conditions:
        parse_condition(c)
    args.outdir.mkdir(parents=True, exist_ok=True)
    image = args.art_image or args.outdir / "art.sif"
    if not image.exists():
        print(f"Downloading ART image -> {image}")
        urllib.request.urlretrieve(ART_IMAGE_URL, image)

    with open(args.assemblies, newline="") as fh:
        fastas = {row["sample"]: Path(row["fasta"]) for row in csv.DictReader(fh)}
    plasmids = plasmid_accessions(args.ground_truth)
    missing = [g for g in args.genomes if g not in fastas or g not in plasmids]
    if missing:
        sys.exit(f"Error: genome(s) not in --assemblies and the ground truth: {missing}")

    workdir = args.outdir / "tmp"
    workdir.mkdir(exist_ok=True)
    rows = []
    for genome in args.genomes:
        for condition in args.conditions:
            sample, r1, r2 = simulate(image, workdir, args.outdir, genome, fastas[genome], plasmids[genome], condition)
            print(f"{sample}: done")
            rows.append((sample, r1.resolve(), r2.resolve()))

    sheet = args.outdir / "sim_samplesheet.csv"
    with open(sheet, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["ID", "R1", "R2", "LongFastQ", "Fast5", "GenomeSize"])
        writer.writerows((s, r1, r2, "NA", "NA", "2.8m") for s, r1, r2 in rows)
    print(f"Wrote {len(rows)} samples to {sheet}")


if __name__ == "__main__":
    main()
