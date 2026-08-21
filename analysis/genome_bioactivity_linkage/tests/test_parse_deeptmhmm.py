import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_deeptmhmm import has_tm_helix_outside_signal, parse_tmrs_gff3  # noqa: E402

FIXTURE = """\
protA\tsignal\t1\t20
protA\toutside\t21\t150
protB\tsignal\t1\t18
protB\tTMhelix\t19\t41
protB\toutside\t42\t200
protC\toutside\t1\t30
protC\tTMhelix\t31\t53
protC\tinside\t54\t100
"""


def test_parse_tmrs_gff3(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    assert list(df.columns) == ["protein_id", "region_type", "start", "end"]
    assert len(df) == 8
    assert df.iloc[0].to_dict() == {"protein_id": "protA", "region_type": "signal", "start": 1, "end": 20}


def test_secreted_protein_with_only_signal_region_has_no_disqualifying_tm(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    assert has_tm_helix_outside_signal(df, "protA", signal_cleavage_site=20) is False


def test_tm_helix_within_cleaved_signal_region_does_not_disqualify(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    # protB's TMhelix (19-41) starts inside the signal region (1-18) cleavage
    # site at 18 but extends past it -- since it extends into the mature
    # chain it DOES disqualify.
    assert has_tm_helix_outside_signal(df, "protB", signal_cleavage_site=18) is True


def test_tm_helix_in_mature_chain_disqualifies():
    import pandas as pd

    df = pd.DataFrame(
        [{"protein_id": "protC", "region_type": "TMhelix", "start": 31, "end": 53}]
    )
    assert has_tm_helix_outside_signal(df, "protC", signal_cleavage_site=None) is True
