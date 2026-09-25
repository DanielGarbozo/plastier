"""Tests for bin/baseline_comparison.py (stage 7d). Run with: pytest tests/python"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "bin" / "baseline_comparison.py"

GENOME = "GCF_000000001.1"
CHROM = "NZ_CP000001.1"
PLASMID = "NZ_CP000002.1"

COLUMNS = ["gene_symbol", "input_sequence_id", "tier", "mobsuite_call", "platon_call", "rfplasmid_call"]
ROWS = [
    # true plasmid, tier High, MOB + RFPlasmid right
    ("g1", PLASMID, "High-confidence plasmid", "plasmid", None, "plasmid"),
    # true plasmid, tier abstains, MOB right, RFPlasmid wrong
    ("g2", PLASMID, "Ambiguous", "plasmid", None, "chromosome"),
    # true chromosome, tier wrongly Moderate, MOB wrong, RFPlasmid has no call
    ("g3", CHROM, "Moderate-confidence plasmid", "plasmid", None, None),
    # true chromosome, everything right, RFPlasmid has no call
    ("g4", CHROM, "Chromosomal", "chromosome", None, None),
]


@pytest.fixture
def compared(tmp_path):
    gt = tmp_path / "gt.csv"
    gt.write_text(f"Genome_Accession,Plasmid_Accession\n{GENOME},{PLASMID}\n")
    pred = tmp_path / f"{GENOME}.tier_resolution.tsv"
    pd.DataFrame(ROWS, columns=COLUMNS).to_csv(pred, sep="\t", index=False)
    out = tmp_path / "cmp.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--predictions", str(pred), "--ground-truth", str(gt), "--output", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return pd.read_csv(out, sep="\t").set_index(["Method", "Class"])


TIER = "Tier framework (High+Moderate = plasmid)"


def test_tier_framework_abstention_costs_recall_and_lowers_coverage(compared):
    plasmid = compared.loc[(TIER, "Plasmid")]
    assert plasmid.Coverage == 0.75
    assert (plasmid.TP, plasmid.FP, plasmid.FN) == (1, 1, 1)  # g2 abstained -> FN
    assert (plasmid.Precision, plasmid.Recall) == (0.5, 0.5)
    chromosome = compared.loc[(TIER, "Chromosome")]
    assert (chromosome.TP, chromosome.FP, chromosome.FN) == (1, 0, 1)  # g3 was called plasmid


def test_single_classifier_is_scored_on_its_own_calls(compared):
    mob = compared.loc[("MOB-suite alone", "Plasmid")]
    assert mob.Coverage == 1.0
    assert (mob.TP, mob.FP, mob.FN) == (2, 1, 0)
    assert round(mob.Precision, 4) == 0.6667 and mob.Recall == 1.0


def test_classifier_with_no_call_counts_as_abstention(compared):
    rfp = compared.loc[("RFPlasmid alone", "Plasmid")]
    assert rfp.Coverage == 0.5  # g3 and g4 have no RFPlasmid call
    assert (rfp.TP, rfp.FP, rfp.FN) == (1, 0, 1)
    rfp_chrom = compared.loc[("RFPlasmid alone", "Chromosome")]
    assert (rfp_chrom.TP, rfp_chrom.FP, rfp_chrom.FN) == (0, 1, 2)


def test_skipped_classifier_is_left_out_not_reported_as_all_abstain(compared):
    methods = {m for m, _ in compared.index}
    assert "Platon alone" not in methods
    assert methods == {TIER, "MOB-suite alone", "RFPlasmid alone"}
