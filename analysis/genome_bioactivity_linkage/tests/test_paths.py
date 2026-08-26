import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import paths  # noqa: E402


def test_species_table_has_both_species():
    assert set(paths.SPECIES) == {"dendrobatidis", "salamandrivorans"}
    bd = paths.SPECIES["dendrobatidis"]
    assert bd["out"] == "Batrachochytrium_dendrobatidis_JEL423"
    assert bd["locustag"] == "FCC698BD"
    bsal = paths.SPECIES["salamandrivorans"]
    assert bsal["out"] == "Batrachochytrium_salamandrivorans_AMFP13"
    assert bsal["locustag"] == "F61BA062"


def test_find_bfd_output_missing_raises_with_expected_path(tmp_path, monkeypatch):
    # Both search roots (BFD's shared run AND this project's local fallback)
    # must be redirected to tmp_path: since 2026-08-25 the real local
    # pfam_hmmscan/... fallback files exist on disk, so a monkeypatched
    # BFD_ROOT alone no longer triggers the missing-path case.
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    monkeypatch.setattr(paths, "GBL_ROOT", tmp_path)
    try:
        paths.find_bfd_output("pfam_hmmscan", "dendrobatidis")
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as e:
        assert "FCC698BD" in str(e)
        assert "pfam_hmmscan" in str(e)


def test_find_bfd_output_finds_locustag_file(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    bucket = tmp_path / "results" / "function" / "pfam_hmmscan" / "00"
    bucket.mkdir(parents=True)
    target = bucket / "FCC698BD.domtblout.gz"
    target.write_bytes(b"")
    found = paths.find_bfd_output("pfam_hmmscan", "dendrobatidis")
    assert found == target


def test_find_bfd_output_with_suffix_disambiguates_domtblout_from_tblout(tmp_path, monkeypatch):
    # Every real pfam_hmmscan locustag bucket contains BOTH a .domtblout.gz
    # AND a .tblout.gz. Without an explicit suffix filter, sorted()[0] picks
    # the domtblout one only by alphabetical luck ('d' < 't') -- construct
    # the ambiguous case explicitly rather than relying on that luck.
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    bucket = tmp_path / "results" / "function" / "pfam_hmmscan" / "00"
    bucket.mkdir(parents=True)
    domtblout = bucket / "FCC698BD.domtblout.gz"
    tblout = bucket / "FCC698BD.tblout.gz"
    domtblout.write_bytes(b"")
    tblout.write_bytes(b"")

    found = paths.find_bfd_output("pfam_hmmscan", "dendrobatidis", suffix=".domtblout.gz")
    assert found == domtblout

    found_tblout = paths.find_bfd_output("pfam_hmmscan", "dendrobatidis", suffix=".tblout.gz")
    assert found_tblout == tblout


def test_find_bfd_output_raises_when_suffix_matches_nothing(tmp_path, monkeypatch):
    # A suffix that matches none of the real files must raise, not silently
    # fall back to the unfiltered match list (which would defeat the whole
    # point of passing suffix -- see test above).
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    bucket = tmp_path / "results" / "function" / "signalp" / "00"
    bucket.mkdir(parents=True)
    (bucket / "FCC698BD.signalp.results.txt.gz").write_bytes(b"")
    try:
        paths.find_bfd_output("signalp", "dendrobatidis", suffix=".signalp.gff3.gz")
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as e:
        assert "FCC698BD" in str(e)
        assert ".signalp.gff3.gz" in str(e)
