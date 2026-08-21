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
