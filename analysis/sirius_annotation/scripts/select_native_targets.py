#!/usr/bin/env python3
"""
Select native SIRIUS targets from the Everything-Bagel aligned feature table.

A target = a feature detected in this project (aligned_features_filled.mgf /
aligned_features.csv) that
  * has MS2 (has_ms2 == True)          -> a block exists in the feature MGF,
  * is singly charged (charge == <--charge>, default 1) -- SIRIUS 6.3.12's
    formula tool cannot process multiply charged precursors, and
  * is NOT already annotated in analysis/sirius_annotation/sirius_annotations.tsv
    (which currently holds the transferred EB annotations).

For the pilot, pass --max-features N to draw a reproducible random subset
(--seed) spread across the feature table; omit it for the full native run.

Output: analysis/sirius_annotation/sirius_native_targets.csv (overwritten),
one row per target with row ID, m/z, RT plus provenance (source_file/source_scan).

Usage:
    python3 scripts/select_native_targets.py [--max-features N] [--seed S] [--charge C]
"""
import argparse
import random
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent.parent
FEATURES = (
    REPO
    / "data"
    / "raw"
    / "gnps2_e9838293_bagel"
    / "nf_output"
    / "feature_finding"
    / "feature_finding_results"
    / "aligned_features.csv"
)
EXISTING = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT_PATH = REPO / "analysis" / "sirius_annotation" / "sirius_native_targets.csv"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-features", type=int, default=None, help="Limit target count (pilot)")
    ap.add_argument("--seed", type=int, default=1234, help="RNG seed for the pilot subset")
    ap.add_argument("--charge", type=int, default=1, help="Only target this charge state (default 1)")
    args = ap.parse_args()

    feat = pd.read_csv(FEATURES, low_memory=False)
    n_all = len(feat)
    print(f"feature table: {n_all} rows", file=sys.stderr)

    already_done = []
    if EXISTING.exists():
        existing = pd.read_csv(EXISTING, sep="\t")
        already_done = sorted(set(existing["row ID"].astype(int)))
    print(f"already annotated: {len(already_done)}", file=sys.stderr)

    cand = feat[
        (feat["has_ms2"] == True)  # noqa: E712
        & (feat["charge"].astype(int) == args.charge)
        & (~feat["row ID"].astype(int).isin(already_done))
    ].copy()
    print(
        f"candidates (has_ms2 & charge={args.charge} & un-annotated): {len(cand)}",
        file=sys.stderr,
    )

    if args.max_features is not None:
        cand = cand.sample(min(args.max_features, len(cand)), random_state=args.seed)
        print(f"pilot subset (seed={args.seed}, max={args.max_features}): {len(cand)}", file=sys.stderr)

    out = cand[["row ID", "row m/z", "row retention time", "source_file", "source_scan"]]
    out = out.drop_duplicates("row ID").sort_values("row ID")
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_PATH, index=False)
    print(f"wrote {len(out)} targets -> {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
