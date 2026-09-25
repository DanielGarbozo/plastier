"""Tests for assets/validation/fetch_reference_genomes.py header rewriting. Run with: pytest tests/python"""

import gzip
import importlib.util
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("fetch_reference_genomes", ROOT / "assets/validation/fetch_reference_genomes.py")
fetch = importlib.util.module_from_spec(spec)
spec.loader.exec_module(fetch)

RAW = b">NZ_CP1.1 Staphylococcus aureus chromosome, complete genome\nACGTAC\nGT\n>NZ_CP2.1 Staphylococcus aureus plasmid p1, complete sequence\nGGGCCC\n"


def test_headers_become_unicycler_style_and_keep_accession_and_sequence():
    out = fetch.unicycler_headers(RAW).decode().splitlines()
    assert out[0] == ">NZ_CP1.1 length=8 depth=1.00x circular=true"
    assert out[3] == ">NZ_CP2.1 length=6 depth=1.00x circular=true"
    assert out[1:3] == ["ACGTAC", "GT"] and out[4] == "GGGCCC"
    assert fetch.sequence_ids(fetch.unicycler_headers(RAW)) == {"NZ_CP1.1", "NZ_CP2.1"}


def test_rewrite_is_idempotent():
    once = fetch.unicycler_headers(RAW)
    assert fetch.unicycler_headers(once) == once


def test_every_replicon_is_marked_circular_so_the_label_is_not_leaked():
    out = fetch.unicycler_headers(RAW).decode()
    assert out.count("circular=true") == 2 == out.count(">")


def test_circularity_coverage_reads_the_rewritten_headers(tmp_path):
    """The point of the rewrite: stage 5c must see circular=True and a depth, not 'no evidence'."""
    fasta = tmp_path / "asm.fna.gz"
    with gzip.open(fasta, "wb") as fh:
        fh.write(fetch.unicycler_headers(RAW))
    sys.path.insert(0, str(ROOT / "bin"))
    import circularity_coverage

    headers = circularity_coverage.load_assembly_headers(str(fasta))
    assert headers == {
        "NZ_CP1.1": {"coverage_ratio": 1.0, "circular": True},
        "NZ_CP2.1": {"coverage_ratio": 1.0, "circular": True},
    }


def test_raw_ncbi_headers_give_no_circularity_evidence(tmp_path):
    """Documents why the rewrite exists: without it the High-confidence tier cannot occur."""
    fasta = tmp_path / "raw.fna"
    fasta.write_bytes(RAW)
    sys.path.insert(0, str(ROOT / "bin"))
    import circularity_coverage

    headers = circularity_coverage.load_assembly_headers(str(fasta))
    assert headers["NZ_CP2.1"] == {"coverage_ratio": None, "circular": False}
