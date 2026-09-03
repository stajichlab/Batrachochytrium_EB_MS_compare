#!/usr/bin/env python3
"""
Build the linked sample-metadata + feature-abundance tables used by every
ordination/pairwise script in this analysis folder.

Source tables:
  - data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/
      feature_finding_results/aligned_features.csv
      (GNPS2 Everything Bagel / FBMN aligned feature table, task
      e983829350de4bb39f278cbf22553247 -- per-row `row ID`/`row m/z`/
      `row retention time` plus one `<sample>.mzML Peak area` column per
      sample; 38,547 features x 123 sample columns)
  - data/metdata/curated_gnps_metadata.tsv (per-sample species/matrix/
    life_stage/use_in_analysis annotation)

Output (analysis/ordination/linked_data/):
  - sample_metadata.csv       one row per FUNGAL analysis sample,
                               condition = f"{matrix}_{life_stage}"
  - feature_abundance.csv.gz  features x samples matrix (row_id, mz, rt,
                               then one column per sample_id), raw peak area
  - blank_metadata.csv        the 30 uninoculated C_liq media blanks, same
                               schema, kept OUT of the analysis matrix
  - blank_abundance.csv.gz    their abundances, for blank-contrast work

This is the Everything-Bagel port of the EB sibling project's MZMINE3-based
build_ordination_table.py; the joined linked_data schema is identical so the
ported ordination/differential scripts run unchanged. Note Bsal filenames use
the plural `_spores` (matching both the curated metadata and the FBMN output);
the `_spore`/`_spores` difference is directional (Bd singular, Bsal plural)
and is handled purely by string-stem matching, never by rewriting.

SAMPLE SELECTION (corrected 2026-09-02 -- this docstring previously asserted
the opposite, and that false claim is why the bug below survived review).

`use_in_analysis == True` does NOT drop the media-blank controls. It drops
IS/QC rows and the B-plate "DO NOT USE" conditioned-media rows only. Of its
90 rows, **30 are uninoculated `*C_liq` media blanks** (`is_C_companion ==
True`). Because blanks exist only for the `liq` matrix, they made every liq
group exactly 50% sterile medium:

    liq_Zoospore   5 fungal + 5 blank      liq_Developed  10 fungal + 10 blank
    spore_*        5/10 fungal + 0 blank

Every liq contrast was therefore a diluted fungal-vs-fungal comparison. The
effect is large and NOT in a consistent direction -- Bd liq Zoospore-vs-
Developed went 536 -> 0 significant when blanks are removed, while Bsal went
54 -> 3,492. We now filter on `is_C_companion != True` and write the blanks
to a separate blank matrix, so they remain available for blank-contrast work
(background_subtraction) without ever entering a biological contrast.

FEATURE SELECTION (added 2026-09-02).

The bagel table's own artifact columns were previously unused, so isotope
peaks, multiply-charged rows, non-default adducts and in-source fragments
were all tested as if they were independent compounds -- inflating feature
counts and padding the BH denominator with non-independent duplicate tests
of the same molecule. Of 38,547 rows: 8,521 (22.1%) are M+1..M+5 isotope
peaks, 3,627 are charge>1, 6,683 are non-default adducts, 471 are ISF. We
keep M+0 & charge==1 & is_default_adduct & !is_isf -> 25,157 features.
Pass --no-artifact-filter to reproduce the pre-2026-09-02 feature universe.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
QUANT = (
    REPO
    / "data" / "raw" / "gnps2_e9838293_bagel"
    / "nf_output" / "feature_finding" / "feature_finding_results"
    / "aligned_features.csv"
)
META = REPO / "data" / "metdata" / "curated_gnps_metadata.tsv"
OUT_DIR = REPO / "analysis" / "ordination" / "linked_data"

AREA_SUFFIX = " Peak area"


def artifact_mask(quant: pd.DataFrame) -> pd.Series:
    """M+0 monoisotopic, singly charged, default adduct, not an in-source fragment."""
    isotope = quant["isotope"].astype("string").fillna("M+0")
    return (
        isotope.isin(["M+0"])
        & (quant["charge"] == 1)
        & (quant["is_default_adduct"] == True)  # noqa: E712
        & (quant["is_isf"] != True)  # noqa: E712
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-artifact-filter", action="store_true",
                    help="keep isotopes/multi-charge/non-default-adduct/ISF rows "
                         "(reproduces the pre-2026-09-02 38,547-feature universe)")
    args = ap.parse_args()

    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta["use_in_analysis"] == True].copy()  # noqa: E712
    # Split the uninoculated C_liq media blanks OUT of the analysis matrix.
    # They are `use_in_analysis == True` but are sterile medium, not samples.
    is_blank = meta["is_C_companion"] == True  # noqa: E712
    blanks = meta[is_blank].copy()
    meta = meta[~is_blank].copy()
    print(f"samples: {len(meta)} fungal, {len(blanks)} media blanks held out",
          file=sys.stderr)
    if len(blanks) != 30:
        sys.exit(f"expected 30 C_liq media blanks in use_in_analysis, found {len(blanks)}")
    for frame in (meta, blanks):
        frame["sample_id"] = frame["filename"].str.replace(".mzML", "", regex=False)
        frame["condition"] = frame["matrix"] + "_" + frame["life_stage"]
        frame["stage_group"] = frame["life_stage"].map(
            {"Zoospore": "Zoospore", "Sporangium": "Developed", "Mature": "Developed"}
        )
        frame["condition_group"] = frame["matrix"] + "_" + frame["stage_group"]

    quant = pd.read_csv(QUANT, low_memory=False)
    n_all = len(quant)
    if not args.no_artifact_filter:
        keep = artifact_mask(quant)
        print(
            f"features: {n_all} -> {int(keep.sum())} after dropping "
            f"{int((quant['isotope'].astype('string').fillna('M+0') != 'M+0').sum())} isotope, "
            f"{int((quant['charge'] != 1).sum())} multi-charge, "
            f"{int((quant['is_default_adduct'] != True).sum())} non-default-adduct, "  # noqa: E712
            f"{int((quant['is_isf'] == True).sum())} ISF rows (overlapping counts)",  # noqa: E712
            file=sys.stderr,
        )
        quant = quant[keep].copy()

    area_cols = {
        c: c[: -len(AREA_SUFFIX)].replace(".mzML", "")
        for c in quant.columns
        if c.endswith(AREA_SUFFIX)
    }
    all_ids = set(meta["sample_id"]) | set(blanks["sample_id"])
    missing = all_ids - set(area_cols.values())
    if missing:
        sys.exit(f"metadata samples missing an area column in feature table: {missing}")

    def build(frame: pd.DataFrame) -> pd.DataFrame:
        cols = [c for c, sid in area_cols.items() if sid in set(frame["sample_id"])]
        out = (
            quant[["row ID", "row m/z", "row retention time"] + cols]
            .rename(columns=area_cols)
            .rename(columns={"row ID": "row_id", "row m/z": "mz",
                             "row retention time": "rt"})
            .fillna(0.0)
        )
        return out

    feat = build(meta)
    blank_feat = build(blanks)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta_cols = [
        "sample_id", "species", "matrix", "life_stage", "condition", "plate",
        "replicate", "is_C_companion", "has_C_companion", "companion_of",
        "timepoint_hrs", "stage_group", "condition_group",
    ]
    meta[meta_cols].to_csv(OUT_DIR / "sample_metadata.csv", index=False)
    feat.to_csv(OUT_DIR / "feature_abundance.csv.gz", index=False)
    blanks[meta_cols].to_csv(OUT_DIR / "blank_metadata.csv", index=False)
    blank_feat.to_csv(OUT_DIR / "blank_abundance.csv.gz", index=False)
    print(
        f"wrote {len(meta)} fungal samples x {len(feat)} features "
        f"(+ {len(blanks)} blanks held out) to {OUT_DIR}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
