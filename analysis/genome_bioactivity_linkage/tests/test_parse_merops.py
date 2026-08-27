import gzip
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_merops import best_merops_hit, load_merops_blasttab  # noqa: E402

# Real MEROPS blastp -outfmt 6 shape: 12 standard tabular BLAST columns, no header.
BLASTTAB = """\
protA\tMER0000002\t34.234\t222\t137\t2\t124\t344\t39\t252\t1.53e-37\t141
protA\tMER0000004\t27.321\t377\t211\t15\t9\t382\t1\t317\t5.33e-29\t112
protB\tMER0208659\t31.010\t287\t182\t7\t162\t438\t1\t281\t8.06e-36\t133
"""

FAMILIES = pd.DataFrame(
    [
        {"mer_id": "MER0000002", "clan": "S01A", "family": "S01.001", "catalytic_type": "serine"},
        {"mer_id": "MER0000004", "clan": "S01A", "family": "S01.152", "catalytic_type": "serine"},
        {"mer_id": "MER0208659", "clan": "MA(E)", "family": "M01.001", "catalytic_type": "metallo"},
    ]
)


def test_load_merops_blasttab_parses_gzipped_outfmt6(tmp_path):
    p = tmp_path / "x.blasttab.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(BLASTTAB)
    df = load_merops_blasttab(p)
    assert set(df["protein_id"]) == {"protA", "protB"}
    assert df["evalue"].dtype.kind == "f"


def test_best_merops_hit_keeps_lowest_evalue_and_annotates_family(tmp_path):
    p = tmp_path / "x.blasttab.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(BLASTTAB)
    df = load_merops_blasttab(p)
    best = best_merops_hit(df, FAMILIES)
    assert len(best) == 2  # one row per protein, not per hit
    row_a = best[best["protein_id"] == "protA"].iloc[0]
    assert row_a["mer_id"] == "MER0000002"  # 1.53e-37 < 5.33e-29
    assert row_a["family"] == "S01.001"
    assert row_a["catalytic_type"] == "serine"
    row_b = best[best["protein_id"] == "protB"].iloc[0]
    assert row_b["catalytic_type"] == "metallo"
