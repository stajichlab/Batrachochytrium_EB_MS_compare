# analysis/genome_bioactivity_linkage/scripts/background_subtraction.py
"""Filter liquid-fraction features by fungal-sample vs C_liq-companion-blank signal.

use_in_analysis == True does NOT exclude the 33 is_C_companion == True
media-blank wells (only IS/QC and B-plate conditioned-media rows are
excluded upstream) -- see Task 2 docstring context in the implementation
plan. This module must be used before any liquid-fraction compound is
treated as a candidate for fungal secretion.
"""
import sys
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

# The smallest positive peak area in the whole bagel table is ~1.4e4, so the
# old _PSEUDOCOUNT = 1.0 was inert: a single spike in one fungal well against
# zeros in the blanks gave log2FC ~ 11 and passed trivially. Use an LOD-scale
# constant instead, so the ratio is damped at the real detection floor.
_PSEUDOCOUNT = 1.4e4

# Minimum number of plate pairs in which a feature must beat its own blank.
_MIN_PAIRED_PLATES = 4


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
    if not blank_samples:
        raise ValueError(
            f"No C_liq companion blank samples found for {species}/{life_stage} in feature table"
        )

    mean_fungal = features[fungal_samples].mean(axis=1)
    mean_blank = features[blank_samples].mean(axis=1)
    log2fc = np.log2((mean_fungal + _PSEUDOCOUNT) / (mean_blank + _PSEUDOCOUNT))

    # PAIRED rule (2026-09-02). The design is 1:1 paired -- A1_liq's blank is
    # A1C_liq, same plate, same replicate, same seed date (`companion_of`) --
    # and the old unpaired group-mean rule threw that away. Measured against a
    # blank-vs-blank null (no fungus on either side, so every pass is false),
    # the mean rule had a ~60% false-pass rate; requiring the feature to beat
    # its OWN plate's blank in >= _MIN_PAIRED_PLATES of 5 plates drops that to
    # ~15% for a quarter of the set size.
    pairs = _plate_pairs(scoped, features)
    if pairs:
        beats = np.zeros(len(features), dtype=int)
        for f_col, b_col in pairs:
            beats += (
                (features[f_col] + _PSEUDOCOUNT) / (features[b_col] + _PSEUDOCOUNT)
                >= min_fc
            ).to_numpy(dtype=int)
        n_pairs = len(pairs)
        passes = beats >= min(_MIN_PAIRED_PLATES, n_pairs)
    else:
        # No resolvable pairing for this stratum -- fall back to the unpaired
        # rule rather than silently passing nothing, and say so.
        print(
            f"[background] {species}/{life_stage}: no plate pairing resolved, "
            "falling back to the unpaired group-mean rule",
            file=sys.stderr,
        )
        beats = np.full(len(features), -1)
        n_pairs = 0
        passes = (log2fc >= np.log2(min_fc)).to_numpy()

    return pd.DataFrame(
        {
            "row_id": features.index,
            "mean_fungal": mean_fungal.values,
            "mean_blank": mean_blank.values,
            "log2fc_fungal_over_blank": log2fc.values,
            "n_plates_beating_blank": beats,
            "n_plate_pairs": n_pairs,
            "passes_background_filter": passes,
        }
    )


def _plate_pairs(scoped: pd.DataFrame, features: pd.DataFrame) -> list[tuple[str, str]]:
    """(fungal_col, blank_col) per plate, matched on (plate, replicate).

    NOT via `companion_of`: on a `*C_liq` blank row that column names the
    *spore* sample of the same well (e.g. `A1C_liq` -> `A1_spore`), not the
    liq sample we want to pair against. Within a species x life_stage x liq
    stratum, (plate, replicate) uniquely identifies the well, so `A1_liq`
    pairs with `A1C_liq` -- same plate, same replicate, same seed date.
    """
    fungal = scoped[~scoped["is_C_companion"]]
    blanks = scoped[scoped["is_C_companion"]]
    by_well = {
        (r["plate"], r["replicate"]): r["filename"] for _, r in blanks.iterrows()
    }
    pairs = []
    for _, r in fungal.iterrows():
        blank_file = by_well.get((r["plate"], r["replicate"]))
        if blank_file and r["filename"] in features.columns and blank_file in features.columns:
            pairs.append((r["filename"], blank_file))
    return pairs
