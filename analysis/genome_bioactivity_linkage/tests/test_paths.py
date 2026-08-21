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
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
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
