"""Tests for assets/validation/plaseval_convert.py. Run with: pytest tests/python"""

import subprocess
import sys
from pathlib import Path

import pandas as pd

SCRIPT = Path(__file__).resolve().parents[2] / "assets" / "validation" / "plaseval_convert.py"

# Column layout and value styles copied from a real MOB-suite contig_report.txt:
# lowercase molecule_type, contig_id = whole FASTA header, "-" for no cluster.
HEADER = "sample_id\tmolecule_type\tprimary_cluster_id\tsecondary_cluster_id\tcontig_id\tsize\tgc\n"


def report(path, sample, rows):
    lines = [f"{sample}\t{mol}\t{cluster}\t-\t{contig}\t{size}\t0.33\n" for mol, cluster, contig, size in rows]
    path.write_text(HEADER + "".join(lines))
    return path


def convert_predicted(tmp_path, *reports):
    out = tmp_path / "pred.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "predicted", "--input", *map(str, reports), "--output", str(out)],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    return pd.read_csv(out, sep="\t", dtype=str)


def test_real_report_style_lowercase_plasmid_and_header_contig_ids(tmp_path):
    r = report(tmp_path / "a.txt", "GCF_1.1", [
        ("chromosome", "-", "NZ_CP1.1 Staphylococcus aureus chromosome", 3000000),
        ("plasmid", "AA849", "NZ_CP2.1 Staphylococcus aureus plasmid 2", 28654),
    ])
    pred = convert_predicted(tmp_path, r)
    assert pred.to_dict("records") == [{"plasmid": "GCF_1.1|AA849", "contig": "NZ_CP2.1", "contig_len": "28654"}]


def test_same_cluster_in_two_genomes_stays_two_bins(tmp_path):
    a = report(tmp_path / "a.txt", "GCF_1.1", [("plasmid", "AA849", "NZ_CP2.1 x", 100)])
    b = report(tmp_path / "b.txt", "GCF_2.1", [("plasmid", "AA849", "NZ_CP9.1 y", 200)])
    pred = convert_predicted(tmp_path, a, b)
    assert sorted(pred["plasmid"]) == ["GCF_1.1|AA849", "GCF_2.1|AA849"]


def test_contigs_in_one_cluster_form_one_bin_and_unclustered_contigs_do_not(tmp_path):
    r = report(tmp_path / "a.txt", "GCF_1.1", [
        ("plasmid", "AB720", "NZ_CP2.1 x", 100),
        ("plasmid", "AB720", "NZ_CP3.1 x", 200),
        ("plasmid", "-", "NZ_CP4.1 x", 300),
        ("plasmid", "-", "NZ_CP5.1 x", 400),
    ])
    pred = convert_predicted(tmp_path, r).set_index("contig")
    assert pred.loc["NZ_CP2.1", "plasmid"] == pred.loc["NZ_CP3.1", "plasmid"] == "GCF_1.1|AB720"
    assert pred.loc["NZ_CP4.1", "plasmid"] == "GCF_1.1|NZ_CP4.1"
    assert pred.loc["NZ_CP5.1", "plasmid"] == "GCF_1.1|NZ_CP5.1"


def test_chromosome_only_report_gives_header_only_output(tmp_path):
    r = report(tmp_path / "a.txt", "GCF_1.1", [("chromosome", "-", "NZ_CP1.1 x", 3000000)])
    pred = convert_predicted(tmp_path, r)
    assert list(pred.columns) == ["plasmid", "contig", "contig_len"] and pred.empty


def test_ground_truth_only_filter(tmp_path):
    gt = tmp_path / "gt.csv"
    gt.write_text(
        "Genome_Accession,Plasmid_Accession,Plasmid_Size_bp\n"
        "GCF_1.1,NZ_CP2.1,100\nGCF_2.1,NZ_CP9.1,200\n"
    )
    out = tmp_path / "gt.tsv"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "ground-truth", "--input", str(gt), "--output", str(out), "--only", "GCF_2.1"],
        capture_output=True, text=True,
    )
    assert proc.returncode == 0, proc.stderr
    assert pd.read_csv(out, sep="\t")["contig"].tolist() == ["NZ_CP9.1"]
