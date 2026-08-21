import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from merge_secretion import predicted_extracellular  # noqa: E402


def test_signal_positive_no_tm_no_gpi_is_extracellular():
    signalp = pd.DataFrame([{"protein_id": "protA", "is_signal_peptide": True, "cleavage_site": 20}])
    deeptmhmm = pd.DataFrame(
        [{"protein_id": "protA", "region_type": "signal", "start": 1, "end": 20}]
    )
    predgpi = pd.DataFrame([{"protein_id": "protA", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protA"].iloc[0]
    assert bool(row["is_extracellular"]) is True


def test_tm_helix_in_mature_chain_excludes():
    signalp = pd.DataFrame([{"protein_id": "protB", "is_signal_peptide": True, "cleavage_site": 18}])
    deeptmhmm = pd.DataFrame(
        [
            {"protein_id": "protB", "region_type": "signal", "start": 1, "end": 18},
            {"protein_id": "protB", "region_type": "TMhelix", "start": 40, "end": 62},
        ]
    )
    predgpi = pd.DataFrame([{"protein_id": "protB", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protB"].iloc[0]
    assert bool(row["is_extracellular"]) is False


def test_gpi_anchor_excludes_even_without_tm():
    signalp = pd.DataFrame([{"protein_id": "protC", "is_signal_peptide": True, "cleavage_site": 22}])
    deeptmhmm = pd.DataFrame(
        [{"protein_id": "protC", "region_type": "signal", "start": 1, "end": 22}]
    )
    predgpi = pd.DataFrame([{"protein_id": "protC", "has_gpi_anchor": True}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protC"].iloc[0]
    assert bool(row["is_extracellular"]) is False


def test_no_signal_peptide_excludes():
    signalp = pd.DataFrame([{"protein_id": "protD", "is_signal_peptide": False, "cleavage_site": None}])
    deeptmhmm = pd.DataFrame(columns=["protein_id", "region_type", "start", "end"])
    predgpi = pd.DataFrame([{"protein_id": "protD", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protD"].iloc[0]
    assert bool(row["is_extracellular"]) is False
