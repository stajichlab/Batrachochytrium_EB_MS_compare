import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from background_subtraction import fungal_over_blank_ratio  # noqa: E402


def _meta_row(filename, species, life_stage, is_companion):
    return {
        "filename": filename,
        "species": species,
        "life_stage": life_stage,
        "is_C_companion": is_companion,
        "matrix": "liq",
        "use_in_analysis": True,
    }


def test_high_fungal_low_blank_passes_filter():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A2_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [1000.0], "A2_liq.mzML": [1200.0], "A1C_liq.mzML": [50.0]},
        index=pd.Index([1], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 1].iloc[0]
    assert row["mean_fungal"] == 1100.0
    assert row["mean_blank"] == 50.0
    assert row["passes_background_filter"] is True or row["passes_background_filter"] == True  # noqa: E712


def test_media_dominated_feature_fails_filter():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [100.0], "A1C_liq.mzML": [90.0]},
        index=pd.Index([2], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 2].iloc[0]
    assert bool(row["passes_background_filter"]) is False


def test_zero_blank_signal_treated_as_pass():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [500.0], "A1C_liq.mzML": [0.0]},
        index=pd.Index([3], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 3].iloc[0]
    assert bool(row["passes_background_filter"]) is True


def test_no_blank_samples_raises_value_error():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [500.0]},
        index=pd.Index([4], name="row_id"),
    )
    with pytest.raises(ValueError):
        fungal_over_blank_ratio(
            features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
        )
