import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_rbh import reciprocal_best_hits  # noqa: E402

# DIAMOND blastp -outfmt 6 default columns:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
FWD = """\
bfd1\tref1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
bfd1\tref2\t80.0\t100\t0\t0\t1\t100\t1\t100\t1e-30\t150
bfd2\tref3\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t180
"""
REV = """\
ref1\tbfd1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
ref3\tbfd2\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t180
ref2\tbfd9\t70.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100
"""


def test_reciprocal_best_hits_only_keeps_mutual_top_hits(tmp_path):
    fwd = tmp_path / "fwd.tsv"
    rev = tmp_path / "rev.tsv"
    fwd.write_text(FWD)
    rev.write_text(REV)
    result = reciprocal_best_hits(fwd, rev)
    pairs = set(zip(result["bfd_protein_id"], result["reference_protein_id"]))
    assert pairs == {("bfd1", "ref1"), ("bfd2", "ref3")}
    assert "bfd1", "ref2" not in pairs  # bfd1's second-best hit, not reciprocal
