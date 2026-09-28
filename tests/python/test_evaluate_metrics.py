"""Tests for bin/evaluate_metrics.py (stage 7b). Run with: pytest tests/python"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "evaluate_metrics.py"
REAL_GROUND_TRUTH = Path(__file__).resolve().parents[2] / "assets" / "validation" / "closed_genome_ground_truth.csv"

GENOME = "GCF_000000001.1"
CHROM = "NZ_CP000001.1"
PLASMID = "NZ_CP000002.1"

GROUND_TRUTH = (
    "Genome_Accession,Plasmid_Accession\n"
    f"{GENOME},{PLASMID}\n"
)


def write_predictions(path, rows):
    pd.DataFrame(rows, columns=["gene_symbol", "input_sequence_id", "tier"]).to_csv(path, sep="\t", index=False)


def run(tmp_path, rows, genome=GENOME, ground_truth=None, extra=()):
    gt = tmp_path / "gt.csv"
    gt.write_text(ground_truth or GROUND_TRUTH)
    pred = tmp_path / f"{genome}.tier_resolution.tsv"
    write_predictions(pred, rows)
    out = tmp_path / "metrics.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--predictions", str(pred), "--ground-truth", str(gt), "--output", str(out), *extra],
        capture_output=True, text=True,
    )
    return proc, out


def metrics(out):
    return pd.read_csv(out, sep="\t").set_index("Tier")


def test_per_tier_counts_and_scores(tmp_path):
    rows = [
        ("blaZ", PLASMID, "High-confidence plasmid"),       # TP
        ("ermC", PLASMID, "High-confidence plasmid"),       # TP
        ("tetM", CHROM, "Moderate-confidence plasmid"),     # FP (actually chromosomal)
        ("mecA", CHROM, "Chromosomal"),                     # TP
        ("aacA", PLASMID, "Chromosomal"),                   # FP (actually plasmid)
        ("qacA", PLASMID, "Ambiguous"),
    ]
    proc, out = run(tmp_path, rows)
    assert proc.returncode == 0, proc.stderr
    m = metrics(out)

    high = m.loc["High-confidence plasmid"]
    assert (high.TP, high.FP, high.FN) == (2, 0, 2)  # 4 true-plasmid ARGs, 2 called High
    assert high.Precision == 1.0 and high.Recall == 0.5

    moderate = m.loc["Moderate-confidence plasmid"]
    assert (moderate.TP, moderate.FP) == (0, 1)
    assert moderate.Precision == 0.0

    chromosomal = m.loc["Chromosomal"]
    assert (chromosomal.TP, chromosomal.FP, chromosomal.FN) == (1, 1, 1)  # 2 true-chromosome ARGs

    ambiguous = m.loc["Ambiguous"]
    assert ambiguous.Total_Calls == 1 and ambiguous.Calls_True_Plasmid == 1
    assert pd.isna(ambiguous.Precision) and pd.isna(ambiguous.Recall) and pd.isna(ambiguous.F1_Score)


def test_same_gene_on_plasmid_and_chromosome_is_scored_per_contig(tmp_path):
    rows = [
        ("qacA", PLASMID, "High-confidence plasmid"),
        ("qacA", CHROM, "Chromosomal"),
    ]
    proc, out = run(tmp_path, rows, extra=["--per-arg-output", str(tmp_path / "per_arg.tsv")])
    assert proc.returncode == 0, proc.stderr
    per_arg = pd.read_csv(tmp_path / "per_arg.tsv", sep="\t").set_index("contig")
    assert per_arg.loc[PLASMID, "true_location"] == "Plasmid"
    assert per_arg.loc[CHROM, "true_location"] == "Chromosome"
    m = metrics(out)
    assert m.loc["High-confidence plasmid", "Precision"] == 1.0
    assert m.loc["Chromosomal", "Precision"] == 1.0


def test_rejects_non_accession_contig_ids(tmp_path):
    proc, _ = run(tmp_path, [("blaZ", "37", "Chromosomal")])
    assert proc.returncode != 0
    assert "not replicon accessions" in proc.stderr


def test_rejects_genome_missing_from_ground_truth(tmp_path):
    proc, _ = run(tmp_path, [("blaZ", CHROM, "Chromosomal")], genome="GCF_999999999.1")
    assert proc.returncode != 0
    assert "not in the ground truth" in proc.stderr


def test_sequence_id_description_is_ignored(tmp_path):
    proc, out = run(tmp_path, [("blaZ", f"{PLASMID} Staphylococcus aureus plasmid, complete sequence", "High-confidence plasmid")])
    assert proc.returncode == 0, proc.stderr
    assert metrics(out).loc["High-confidence plasmid", "TP"] == 1


def test_real_ground_truth_schema_loads(tmp_path):
    """Regression: the script once read columns that do not exist in the real ground-truth CSV."""
    gt = pd.read_csv(REAL_GROUND_TRUTH, dtype=str)
    genome, plasmid = gt.iloc[0][["Genome_Accession", "Plasmid_Accession"]]
    proc, out = run(
        tmp_path,
        [("blaZ", plasmid, "High-confidence plasmid")],
        genome=genome,
        ground_truth=REAL_GROUND_TRUTH.read_text(),
    )
    assert proc.returncode == 0, proc.stderr
    assert metrics(out).loc["High-confidence plasmid", "TP"] == 1


def test_genome_id_flag_requires_single_file(tmp_path):
    proc, _ = run(tmp_path, [("blaZ", CHROM, "Chromosomal")], extra=["--genome-id", GENOME, "--predictions", "a", "b"])
    assert proc.returncode != 0
