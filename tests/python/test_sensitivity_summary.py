"""Tests for bin/sensitivity_summary.py (stage 7e, diagnostic only). Run with: pytest tests/python"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "sensitivity_summary.py"
HIGH, MOD, CHROM, AMBIG = "High-confidence plasmid", "Moderate-confidence plasmid", "Chromosomal", "Ambiguous"


def tier_file(path, rows):
    pd.DataFrame(rows, columns=["gene_symbol", "tier"]).to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def result(tmp_path):
    ref = tier_file(tmp_path / "GCF_1.1.tier_resolution.tsv", [
        ("blaZ", HIGH), ("ermC", CHROM), ("mecA", CHROM), ("tetM", MOD), ("qacA", CHROM), ("qacA", HIGH),
    ])
    sim = tier_file(tmp_path / "GCF_1.1_d20_q0_p1.tier_resolution.tsv", [
        ("blaZ", HIGH),   # same
        ("ermC", AMBIG),  # chromosomal -> ambiguous
        # mecA: lost
        ("tetM", CHROM),  # plasmid -> chromosome: flipped
        ("qacA", CHROM),  # gene in two tiers before, one now
    ])
    empty = tier_file(tmp_path / "GCF_1.1_d5_q-15_p1.tier_resolution.tsv", [])  # a run that found nothing
    out, detail = tmp_path / "s.tsv", tmp_path / "d.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--reference", str(ref), "--simulated", str(sim), str(empty),
         "--output", str(out), "--detail-output", str(detail)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return pd.read_csv(out, sep="\t"), pd.read_csv(detail, sep="\t")


def outcomes(detail, sample):
    return detail[detail["sample"] == sample].set_index("gene_symbol")["Outcome"].to_dict()


def test_each_kind_of_change_is_classified(result):
    _, detail = result
    assert outcomes(detail, "GCF_1.1_d20_q0_p1") == {
        "blaZ": "same", "ermC": "to_ambiguous", "mecA": "lost", "tetM": "flipped", "qacA": "changed_other",
    }


def test_run_that_found_nothing_counts_every_gene_as_lost(result):
    _, detail = result
    assert set(outcomes(detail, "GCF_1.1_d5_q-15_p1").values()) == {"lost"}
    assert len(outcomes(detail, "GCF_1.1_d5_q-15_p1")) == 5


def test_summary_is_per_condition_and_reference_tier(result):
    summary, _ = result
    row = summary[(summary.Depth == 20) & (summary.Reference_Tier == CHROM)].iloc[0]
    assert (row.N_Genes, row.Recovered, row.Same_Tier, row.To_Ambiguous) == (2, 1, 0, 1)  # ermC, mecA
    assert row.Recovered_Pct == 50.0
    mixed = summary[summary.Reference_Tier.str.startswith("Mixed")]
    assert set(mixed.N_Genes) == {1}


def test_every_summary_row_is_labelled_diagnostic(result):
    summary, _ = result
    assert summary["Note"].str.contains("DIAGNOSTIC ONLY").all()


def test_unrecognised_sample_name_is_rejected(tmp_path):
    ref = tier_file(tmp_path / "GCF_1.1.tier_resolution.tsv", [("blaZ", HIGH)])
    bad = tier_file(tmp_path / "not_a_condition.tier_resolution.tsv", [("blaZ", HIGH)])
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--reference", str(ref), "--simulated", str(bad), "--output", str(tmp_path / "o.tsv")],
        capture_output=True, text=True,
    )
    assert proc.returncode != 0 and "not named" in proc.stderr


def test_assembler_suffix_on_sample_id_is_accepted(tmp_path):
    ref = tier_file(tmp_path / "GCF_1.1.tier_resolution.tsv", [("blaZ", HIGH)])
    sim = tier_file(tmp_path / "GCF_1.1_d20_q0_p1-unicycler.tier_resolution.tsv", [("blaZ", HIGH)])
    out = tmp_path / "o.tsv"
    proc = subprocess.run([sys.executable, str(SCRIPT), "--reference", str(ref), "--simulated", str(sim), "--output", str(out)],
                          capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    assert pd.read_csv(out, sep="\t").iloc[0].Same_Tier == 1


def test_sample_with_no_tier_file_counts_as_every_gene_lost(tmp_path):
    """An assembly that fails writes nothing; it must stay in the denominator."""
    ref = tier_file(tmp_path / "GCF_1.1.tier_resolution.tsv", [("blaZ", HIGH), ("ermC", CHROM)])
    ok = tier_file(tmp_path / "GCF_1.1_d50_q0_p1-unicycler.tier_resolution.tsv", [("blaZ", HIGH), ("ermC", CHROM)])
    sheet = tmp_path / "sheet.csv"
    sheet.write_text("ID,R1,R2\nGCF_1.1_d50_q0_p1,a,b\nGCF_1.1_d20_q-15_p1,a,b\n")  # second one never assembled
    out, detail = tmp_path / "o.tsv", tmp_path / "d.tsv"
    proc = subprocess.run([sys.executable, str(SCRIPT), "--reference", str(ref), "--simulated", str(ok), "--samplesheet", str(sheet),
                           "--output", str(out), "--detail-output", str(detail)], capture_output=True, text=True)
    assert proc.returncode == 0, proc.stderr
    d = pd.read_csv(detail, sep="\t")
    failed = d[d["Quality_Shift"] == -15]
    assert len(failed) == 2 and set(failed.Outcome) == {"lost"}
    assert set(d[d["Quality_Shift"] == 0].Outcome) == {"same"}
