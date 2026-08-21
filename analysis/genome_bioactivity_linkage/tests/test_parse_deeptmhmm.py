import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_deeptmhmm import has_tm_helix_outside_signal, parse_tmrs_gff3  # noqa: E402

# Real-shaped excerpt trimmed from a live DeepTMHMM run's TMRs.gff3
# (/rhome/jstajich/projects/nf/nf_funannotate1/tests/output/deeptmhmm_gpu_test/
# 27596677/TMRs.gff3): data lines tab-padded to 8 fields, records separated by
# bare "//" lines (not "#" comments). B2ANF9.1 is copied verbatim from that
# file; protD is a constructed (but real-shaped) record added to exercise a
# TMhelix that falls entirely within the signal-peptide region.
FIXTURE = (
    "##gff-version 3\n"
    "# B2ANF9.1 Number of predicted TMRs: 4\n"
    "B2ANF9.1\tinside\t1\t36\t\t\t\t\n"
    "B2ANF9.1\tTMhelix\t37\t57\t\t\t\t\n"
    "B2ANF9.1\toutside\t58\t78\t\t\t\t\n"
    "B2ANF9.1\tTMhelix\t79\t97\t\t\t\t\n"
    "B2ANF9.1\tinside\t98\t105\t\t\t\t\n"
    "//\n"
    "# P40578.1 Number of predicted TMRs: 1\n"
    "P40578.1\tinside\t1\t1036\t\t\t\t\n"
    "P40578.1\tTMhelix\t1037\t1055\t\t\t\t\n"
    "P40578.1\toutside\t1056\t1113\t\t\t\t\n"
    "//\n"
    "# protD.1 Number of predicted TMRs: 1 (constructed: TMhelix confined\n"
    "# entirely within the cleaved signal-peptide region)\n"
    "protD\tsignal\t1\t25\t\t\t\t\n"
    "protD\tTMhelix\t10\t20\t\t\t\t\n"
    "protD\toutside\t26\t200\t\t\t\t\n"
    "//\n"
)


def test_parse_tmrs_gff3(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    assert list(df.columns) == ["protein_id", "region_type", "start", "end"]
    # "##gff-version 3", "# ..." comment lines, and "//" separators must all
    # be skipped; only genuine data lines counted (5 + 3 + 3 = 11).
    assert len(df) == 11
    assert df.iloc[0].to_dict() == {"protein_id": "B2ANF9.1", "region_type": "inside", "start": 1, "end": 36}


def test_tm_helix_in_mature_chain_disqualifies(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    # B2ANF9.1 has TMhelix rows extending well past any plausible cleavage
    # site -- disqualifying.
    assert has_tm_helix_outside_signal(df, "B2ANF9.1", signal_cleavage_site=36) is True


def test_tm_helix_confined_within_signal_region_does_not_disqualify(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    # protD's TMhelix (10-20) falls entirely within the cleaved signal
    # region (cleavage_site=25) -- end <= cleavage_site -- so it must NOT
    # disqualify the protein from being predicted extracellular. This is the
    # true positive-non-disqualifying overlap case: unlike the "no TMhelix
    # rows at all" fixture, this protein genuinely HAS a TMhelix row, just
    # one that is fully upstream of the cleavage site.
    assert has_tm_helix_outside_signal(df, "protD", signal_cleavage_site=25) is False


def test_secreted_protein_with_only_signal_region_has_no_disqualifying_tm():
    import pandas as pd

    df = pd.DataFrame(
        [{"protein_id": "protA", "region_type": "signal", "start": 1, "end": 20}]
    )
    assert has_tm_helix_outside_signal(df, "protA", signal_cleavage_site=20) is False


def test_no_tm_rows_and_no_cleavage_site_is_not_disqualifying():
    import pandas as pd

    df = pd.DataFrame(columns=["protein_id", "region_type", "start", "end"])
    assert has_tm_helix_outside_signal(df, "protX", signal_cleavage_site=None) is False
