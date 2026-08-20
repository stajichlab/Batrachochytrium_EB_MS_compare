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
  - sample_metadata.csv       one row per analysis sample (use_in_analysis
                               == True), condition = f"{matrix}_{life_stage}"
  - feature_abundance.csv.gz  features x samples matrix (row_id, mz, rt,
                               then one column per sample_id), raw peak area

This is the Everything-Bagel port of the EB sibling project's MZMINE3-based
build_ordination_table.py; the joined linked_data schema is identical so the
ported ordination/differential scripts run unchanged. Note Bsal filenames use
the plural `_spores` (matching both the curated metadata and the FBMN output);
the `_spore`/`_spores` difference is directional (Bd singular, Bsal plural)
and is handled purely by string-stem matching, never by rewriting.

use_in_analysis == True already drops IS/QC rows and media-blank controls
(see data/metdata/curated_gnps_metadata.tsv provenance); "not applicable"
species rows (QC/IS) are excluded here as well since they carry no
species/state label to ordinate by.
"""
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


def main():
    meta = pd.read_csv(META, sep="\t")
    meta = meta[meta["use_in_analysis"] == True].copy()  # noqa: E712
    meta["sample_id"] = meta["filename"].str.replace(".mzML", "", regex=False)
    meta["condition"] = meta["matrix"] + "_" + meta["life_stage"]
    meta["stage_group"] = meta["life_stage"].map(
        {"Zoospore": "Zoospore", "Sporangium": "Developed", "Mature": "Developed"}
    )
    meta["condition_group"] = meta["matrix"] + "_" + meta["stage_group"]

    quant = pd.read_csv(QUANT, low_memory=False)
    area_cols = {
        c: c[: -len(AREA_SUFFIX)].replace(".mzML", "")
        for c in quant.columns
        if c.endswith(AREA_SUFFIX)
    }
    keep_samples = [c for c, sid in area_cols.items() if sid in set(meta["sample_id"])]
    missing = set(meta["sample_id"]) - set(area_cols.values())
    if missing:
        sys.exit(f"metadata samples missing an area column in feature table: {missing}")

    feat = (
        quant[["row ID", "row m/z", "row retention time"] + keep_samples]
        .rename(columns=area_cols)
    )
    feat = feat.rename(
        columns={
            "row ID": "row_id",
            "row m/z": "mz",
            "row retention time": "rt",
        }
    ).fillna(0.0)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    meta_cols = [
        "sample_id", "species", "matrix", "life_stage", "condition", "plate",
        "replicate", "is_C_companion", "has_C_companion", "companion_of",
        "timepoint_hrs", "stage_group", "condition_group",
    ]
    meta[meta_cols].to_csv(OUT_DIR / "sample_metadata.csv", index=False)
    feat.to_csv(OUT_DIR / "feature_abundance.csv.gz", index=False)
    print(
        f"wrote {len(meta)} samples x {len(feat)} features to {OUT_DIR}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
