# analysis/genome_bioactivity_linkage/scripts/background_subtraction.py
"""Filter liquid-fraction features by fungal-sample vs C_liq-companion-blank signal.

use_in_analysis == True does NOT exclude the 33 is_C_companion == True
media-blank wells (only IS/QC and B-plate conditioned-media rows are
excluded upstream) -- see Task 2 docstring context in the implementation
plan. This module must be used before any liquid-fraction compound is
treated as a candidate for fungal secretion.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from paths import REPO_ROOT

METADATA_PATH = REPO_ROOT / "data" / "metdata" / "curated_gnps_metadata.tsv"
FEATURES_PATH = (
    REPO_ROOT
    / "data"
    / "raw"
    / "gnps2_e9838293_bagel"
    / "nf_output"
    / "feature_finding"
    / "feature_finding_results"
    / "aligned_features.csv"
)

_PSEUDOCOUNT = 1.0


def load_metadata() -> pd.DataFrame:
    return pd.read_csv(METADATA_PATH, sep="\t")


def load_feature_intensities() -> pd.DataFrame:
    df = pd.read_csv(FEATURES_PATH)
    id_col = "row ID" if "row ID" in df.columns else "row_id"
    df = df.set_index(id_col)
    df.index.name = "row_id"
    peak_area_cols = {c: c.replace(" Peak area", "") for c in df.columns if c.endswith(" Peak area")}
    return df[list(peak_area_cols)].rename(columns=peak_area_cols)


def fungal_over_blank_ratio(
    features: pd.DataFrame,
    meta: pd.DataFrame,
    species: str,
    life_stage: str,
    min_fc: float = 2.0,
) -> pd.DataFrame:
    scoped = meta[
        (meta["species"] == species)
        & (meta["life_stage"] == life_stage)
        & (meta["matrix"] == "liq")
        & (meta["use_in_analysis"] == True)  # noqa: E712
    ]
    fungal_samples = [f for f in scoped.loc[~scoped["is_C_companion"], "filename"] if f in features.columns]
    blank_samples = [f for f in scoped.loc[scoped["is_C_companion"], "filename"] if f in features.columns]
    if not fungal_samples:
        raise ValueError(f"No fungal liq samples found for {species}/{life_stage} in feature table")

    mean_fungal = features[fungal_samples].mean(axis=1)
    mean_blank = features[blank_samples].mean(axis=1) if blank_samples else pd.Series(0.0, index=features.index)

    log2fc = np.log2((mean_fungal + _PSEUDOCOUNT) / (mean_blank + _PSEUDOCOUNT))
    passes = log2fc >= np.log2(min_fc)

    return pd.DataFrame(
        {
            "row_id": features.index,
            "mean_fungal": mean_fungal.values,
            "mean_blank": mean_blank.values,
            "log2fc_fungal_over_blank": log2fc.values,
            "passes_background_filter": passes.values,
        }
    )
