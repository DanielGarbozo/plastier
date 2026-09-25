"""
Regression test for stage 5's per-hit join (bin/*.py signals -> bin/tier_resolution.py).

Found on real data (GCF_045345275.1): Tn554 (ant(9)-Ia + erm(A)) sits twice on the
chromosome, once inside an SCCmec cassette and once outside. The signal tables were
keyed on (gene, contig) only, so each such gene came out 2**5 = 32 times and one
copy's SCCmec verdict was paired with the other copy's. Run with: pytest tests/python
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

BIN = Path(__file__).resolve().parents[2] / "bin"
SAMPLE = "S1"
CONTIG = "chrA"
CONTIG_LEN = 2000  # above the 1000 bp default length floor


def run(script, *args):
    proc = subprocess.run([sys.executable, str(BIN / script), *map(str, args)], capture_output=True, text=True)
    assert proc.returncode == 0, f"{script} failed:\n{proc.stderr}"


def build_inputs(tmp_path, hits):
    """hits: [(gene, start, stop)] - all on CONTIG. One SCCmec element spans 5-120."""
    ham = tmp_path / "ham.tsv"
    pd.DataFrame(
        [(SAMPLE, CONTIG, gene, start, stop) for gene, start, stop in hits],
        columns=["input_file_name", "input_sequence_id", "gene_symbol", "input_gene_start", "input_gene_stop"],
    ).to_csv(ham, sep="\t", index=False)

    mob = tmp_path / "contig_report.txt"
    mob.write_text(
        "sample_id\tmolecule_type\tcontig_id\trep_type(s)\trelaxase_type(s)\tpredicted_mobility\n"
        f"{SAMPLE}\tchromosome\t{CONTIG}\t-\t-\t-\n"
    )
    fasta = tmp_path / "asm.fa"
    fasta.write_text(f">{CONTIG} length={CONTIG_LEN} depth=1.0x\n{'ACGT' * (CONTIG_LEN // 4)}\n")
    sccmec = tmp_path / "sccmec.tsv"
    sccmec.write_text(
        "Input_File\tStatus\tContig\tAttL_Start\tAttL_End\tAttR_Start\tAttR_End\n"
        f"{SAMPLE}\textracted\t{CONTIG}\t5\t20\t100\t120\n"
    )
    return ham, mob, fasta, sccmec


def resolve(tmp_path, hits):
    ham, mob, fasta, sccmec = build_inputs(tmp_path, hits)
    common = ["--hamronization", ham, "--sample-id", SAMPLE]
    out = {name: tmp_path / f"{name}.tsv" for name in ("agree", "marker", "cov", "length", "sccmec", "tier")}
    run("classifier_agreement.py", *common, "--mobsuite", mob, "--output", out["agree"])
    run("mobility_markers.py", *common, "--mobsuite", mob, "--output", out["marker"])
    run("circularity_coverage.py", *common, "--assembly", fasta, "--output", out["cov"])
    run("contig_length_floor.py", *common, "--assembly", fasta, "--output", out["length"])
    run("sccmec_override.py", *common, "--sccmecextractor", sccmec, "--output", out["sccmec"])
    run(
        "tier_resolution.py",
        "--classifier-agreement", out["agree"], "--mobility-markers", out["marker"],
        "--circularity-coverage", out["cov"], "--contig-length-floor", out["length"],
        "--sccmec-override", out["sccmec"], "--output", out["tier"],
    )
    return pd.read_csv(out["tier"], sep="\t", dtype=str)


def test_repeated_gene_is_resolved_per_copy(tmp_path):
    tier = resolve(tmp_path, [("erm(A)", 10, 100), ("erm(A)", 400, 490)])

    assert len(tier) == 2, f"expected one row per copy, got {len(tier)}:\n{tier}"
    by_start = tier.set_index("input_gene_start")
    assert by_start.loc["10", "sccmec_override"] == "True"     # inside the cassette
    assert by_start.loc["400", "sccmec_override"] == "False"   # outside it
    assert by_start.loc["10", "sccmec_element_span"] == "5-120"


def test_same_hit_reported_by_several_tools_collapses_to_one_row(tmp_path):
    tier = resolve(tmp_path, [("blaZ", 200, 300)] * 3)
    assert len(tier) == 1


def test_no_row_multiplication_with_many_repeats(tmp_path):
    hits = [("erm(A)", 10 + 50 * i, 40 + 50 * i) for i in range(4)] + [("blaZ", 500, 590)]
    tier = resolve(tmp_path, hits)
    assert len(tier) == 5


def test_tier_output_keeps_leading_columns_and_adds_coordinates(tmp_path):
    tier = resolve(tmp_path, [("blaZ", 200, 300)])
    assert list(tier.columns[:4]) == ["gene_symbol", "input_sequence_id", "input_gene_start", "input_gene_stop"]
    assert tier.loc[0, "tier"] == "Chromosomal"
