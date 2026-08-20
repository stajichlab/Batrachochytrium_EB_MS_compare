#!/usr/bin/env python3
"""
Transfer SIRIUS annotations from the sibling EB project's run onto THIS
project's (Everything-Bagel / GNPS2 task e983829350de4bb39f278cbf22553247)
feature table, accumulated into analysis/sirius_annotation/sirius_annotations.tsv
keyed by THIS project's `row ID` (feature id).

Why a transfer is needed (and exact, not approximate):
  The two projects process the SAME MassIVE deposit (MSV000090464) with
  DIFFERENT feature finders:
    ../EB  -> MZMine3 `gnps_ms2.mgf`, 4,107 features, SIRIUS keyed to those row IDs.
    this   -> GNPS2 Everything Bagel, 38,547 features (`aligned_features.csv`,
              SCANS = feature id in `aligned_features_filled.mgf`).
  Feature IDs therefore do not correspond, but the underlying LC-MS signal is
  identical: true counterpart features sit at ~0 ppm precursor m/z with RT
  within a few tenths of a minute. So the join is m/z + RT candidate search
  followed by MS2 spectral-cosine disambiguation against both local MGFs.

Decision tree (per EB SIRIUS-annotated feature, tested at defaults):
  1. candidates = current features with |ppm(mz)| <= --ppm AND |dRT| <= --rt
  2. score each scoring candidate by MS2 cosine (fragment tol 0.05 Da,
     --frag-tol) against the EB consensus spectrum; candidates whose current
     feature has no positive-peak MS2 (gap-filled placeholders) cannot be
     scored.
  3. no feature in the m/z+RT window    -> status no_match          (unresolved)
     candidates exist but none scoreable-> status no_ms2_candidate (unresolved)
     best cosine <  --cos-min           -> status low_cosine       (unresolved)
     best cosine >= --cos-min           -> ASSIGN to best candidate:
         single scoreable           -> match_class unique
         2nd scoreable also >= 0.7  -> match_class ms2_tie  (near-duplicate
                                       features of the same compound; pick the
                                       highest-cosine one, keep provenance)
         otherwise                  -> match_class ms2_winner
  Defaults (PPM=10, RT=0.5, COS_MIN=0.7, FRAG_TOL=0.05) assign
  ~76% of the 2,860 EB-annotated features; the residual ~24% (no_match,
  no_ms2_candidate, low_cosine) stay unresolved and are listed in
  sirius_transfer_map.tsv for a later manual pass.

Merging the next (native) SIRIUS round:
  When SIRIUS has been re-run on THIS project's own `aligned_features_filled.mgf`
  (SCANS == feature id, same 7-table merged layout as EB's bundle), import it
  with:
    python3 analysis/sirius_annotation/scripts/import_sirius_transfer.py \
      --native-merged <merged_dir> [--native-label <label>]
  Native rows (annotation_origin=native, keyed directly to this project's
  feature ids -- no mapping) are added to the same accumulated table and WON
  over transferred rows for the same feature id (precedence: native >
  transferred; then structure-hit > none; then higher structure confidence),
  exactly mirroring ../EB/scripts/import_sirius_annotations.py. Transferred
  rows survive only for feature ids the native run did not annotate (e.g. the
  multiply-charged ~886 features SIRIUS skips, and gap-filled placeholders).
  Every contributing source run is preserved in the provenance columns.

Usage:
  # default: transfer from ../EB (as configured below)
  python3 analysis/sirius_annotation/scripts/import_sirius_transfer.py
  # later: add this project's own native run and re-merge
  python3 analysis/sirius_annotation/scripts/import_sirius_transfer.py \
      --native-merged path/to/native/merged --native-label native-e9838293-bagel
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent.parent   # <- scripts/..4 = repo root
ANALYSIS = REPO / "analysis" / "sirius_annotation"
SCRIPTS = ANALYSIS / "scripts"

# ---- inputs (sibling EB project, consumed read-only) ----
EB = Path("/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB")
EB_SIRIUS_ANNOT = EB / "analysis/sirius_annotation/sirius_annotations.tsv"
EB_TARGETS_CSV = EB / "analysis/sirius_annotation/sirius_targets.csv"
EB_TARGETS_MGF = EB / "analysis/sirius_annotation/sirius_targets.mgf"
EB_SOURCE_RUN = "EB-MSV000090464-gnps_ms2"

# ---- inputs (this project's Everything-Bagel bundle, immutable) ----
FEAT_CSV = (
    REPO / "data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/"
    "feature_finding_results/aligned_features.csv"
)
FEAT_MGF = (
    REPO / "data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/"
    "aligned_features_filled.mgf"
)

# ---- outputs ----
OUT_ANNOT = ANALYSIS / "sirius_annotations.tsv"
OUT_MAP = ANALYSIS / "sirius_transfer_map.tsv"

# SIRIUS chemistry columns carried over from the EB run (same vocabulary as
# ../EB/scripts/import_sirius_annotations.py)
CHEM_COLS = [
    "sirius_formula", "sirius_adduct", "sirius_structure_name",
    "sirius_structure_smiles", "sirius_structure_confidence",
    "sirius_npc_pathway", "sirius_npc_class", "sirius_classyfire_class",
]

# Provenance / match-tracing columns attached to every accumulated row.
PROV_COLS = [
    "sirius_source_feature_id",   # feature id inside the source SIRIUS run (EB row id, or native feature id)
    "sirius_source_run",          # label of the SIRIUS run the prediction came from
    "annotation_origin",          # 'transferred' (mapped from EB) | 'native' (direct)
    "match_status",               # assigned | no_match | no_ms2_candidate | low_cosine | native
    "match_class",                # unique | ms2_winner | ms2_tie | native
    "n_candidates",               # current features within m/z+RT window
    "source_mz", "source_rt",     # EB (source-run) precursor m/z, RT (min)
    "feature_mz", "feature_rt",   # chosen this-project feature m/z, RT (min)
    "ppm_error", "rt_delta_min",  # source run -> chosen feature
    "ms2_cosine",                 # MS2 spectral cosine of the chosen assignment
    "n_sirius_hits",              # distinct source-run feature ids landing on this feature
    "sirius_hit_ids",             # semicolon-joined source-run feature ids
    "sirius_hit_formulas",        # semicolon-joined formulas carried by those hits
    "merged_conflict",            # True when the collapsed hits disagree on formula
]


# --------------------------------------------------------------------------- #
# MGF + spectral helpers (stdlib/numpy only -- no pyteomics dependency)
# --------------------------------------------------------------------------- #
def parse_mgf(path: Path):
    """SCANS -> (m/z array, intensity array) for spectra that carry at least
    one positive-intensity peak. Empty / gap-filled placeholder spectra are
    omitted (their SCANS are simply absent from the returned dict)."""
    specs: dict[int, tuple[np.ndarray, np.ndarray]] = {}
    cur = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line == "BEGIN IONS":
                cur = {"s": None, "peaks": []}
            elif line == "END IONS":
                if cur is not None and cur["s"] is not None and cur["peaks"]:
                    mz = np.array([p[0] for p in cur["peaks"]])
                    it = np.array([p[1] for p in cur["peaks"]], dtype=float)
                    keep = (it > 0) & np.isfinite(it)
                    if keep.any():
                        specs[cur["s"]] = (mz[keep], it[keep])
                cur = None
            elif "=" in line and cur is not None:
                k, v = line.split("=", 1)
                if k.strip() == "SCANS":
                    cur["s"] = int(v.strip())
            elif cur is not None:
                p = line.split()
                if len(p) >= 2:
                    try:
                        cur["peaks"].append((float(p[0]), float(p[1])))
                    except ValueError:
                        pass
    return specs


def ms2_cosine(a, b, tol=0.05) -> float:
    """Cosine of two peak lists, peaks matched when |mz1-mz2| <= tol,
    intensities L2-normalised before the dot product (plain spectral cosine).
    Returns 0.0 for any empty/absent spectrum."""
    mza, ia = a
    mzb, ib = b
    if len(mza) == 0 or len(mzb) == 0:
        return 0.0
    ia = ia / np.sqrt(np.sum(ia ** 2))
    ib = ib / np.sqrt(np.sum(ib ** 2))
    oa = np.argsort(mza)
    ob = np.argsort(mzb)
    mza, ia = mza[oa], ia[oa]
    mzb, ib = mzb[ob], ib[ob]
    i = j = 0
    dot = 0.0
    while i < len(mza) and j < len(mzb):
        d = mza[i] - mzb[j]
        if d < -tol:
            i += 1
        elif d > tol:
            j += 1
        else:
            dot += ia[i] * ib[j]
            i += 1
            j += 1
    return float(dot)


# --------------------------------------------------------------------------- #
# Native SIRIUS import (mirrors ../EB/scripts/import_sirius_annotations.py)
# --------------------------------------------------------------------------- #
def load_merged(merged_dir: Path, fname: str) -> pd.DataFrame | None:
    path = merged_dir / fname
    if not path.exists():
        print(f"skip: {path} not found", file=sys.stderr)
        return None
    df = pd.read_csv(path, sep="\t")
    df = df.rename(columns={"mappingFeatureId": "row ID"})
    df["row ID"] = df["row ID"].astype(int)
    return df


def distill_native_run(merged_dir: Path, label: str) -> pd.DataFrame:
    """One row per feature id from a finished/merged native SIRIUS run --
    top formula, top structure, CANOPUS class -- keyed on `row ID` (== the
    feature id of THIS project's own MGF)."""
    formula = load_merged(merged_dir, "formula_identifications.tsv")
    if formula is None:
        sys.exit("no formula_identifications.tsv in --native-merged dir "
                 "-- has SIRIUS finished/been merged?")
    structure = load_merged(merged_dir, "structure_identifications.tsv")
    canopus = load_merged(merged_dir, "canopus_structure_summary.tsv")
    if canopus is None:
        canopus = load_merged(merged_dir, "canopus_formula_summary.tsv")

    best_formula = (
        formula[formula["formulaRank"] == 1][["row ID", "molecularFormula", "adduct"]]
        .drop_duplicates("row ID")
        .rename(columns={"molecularFormula": "sirius_formula", "adduct": "sirius_adduct"})
    )
    parts = [best_formula]

    if structure is not None:
        best_structure = (
            structure[structure["structurePerIdRank"] == 1]
            [["row ID", "name", "smiles", "ConfidenceScoreExact"]]
            .drop_duplicates("row ID")
            .rename(columns={"name": "sirius_structure_name",
                             "smiles": "sirius_structure_smiles",
                             "ConfidenceScoreExact": "sirius_structure_confidence"})
        )
        best_structure["sirius_structure_confidence"] = (
            best_structure["sirius_structure_confidence"].replace(
                [float("inf"), float("-inf")], pd.NA))
        parts.append(best_structure)

    if canopus is not None:
        map_cols = {"row ID": "row ID",
                    "NPC#pathway": "sirius_npc_pathway",
                    "NPC#class": "sirius_npc_class",
                    "ClassyFire#class": "sirius_classyfire_class"}
        avail = [c for c in map_cols if c in canopus.columns]
        parts.append(canopus[avail].drop_duplicates("row ID").rename(columns=map_cols))

    out = parts[0]
    for p in parts[1:]:
        out = out.merge(p, on="row ID", how="outer")

    # provenance columns for native rows (no mapping involved)
    out["sirius_source_feature_id"] = out["row ID"]
    out["sirius_source_run"] = label
    out["annotation_origin"] = "native"
    out["match_status"] = "native"
    out["match_class"] = "native"
    out["n_candidates"] = 1
    out["source_mz"] = np.nan
    out["source_rt"] = np.nan
    out["feature_mz"] = np.nan
    out["feature_rt"] = np.nan
    out["ppm_error"] = np.nan
    out["rt_delta_min"] = np.nan
    out["ms2_cosine"] = np.nan
    return out


# --------------------------------------------------------------------------- #
# Transfer core
# --------------------------------------------------------------------------- #
def load_current_features():
    """(feature_id, mz, rt, has_ms2) arrays from aligned_features.csv."""
    import csv as _csv
    fids, mz, rt, ms2 = [], [], [], []
    with open(FEAT_CSV) as f:
        for row in _csv.DictReader(f):
            fids.append(int(row["row ID"]))
            mz.append(float(row["row m/z"]))
            rt.append(float(row["row retention time"]))
            ms2.append(row["has_ms2"] == "true")
    return (np.array(fids), np.array(mz), np.array(rt), np.array(ms2))


def load_eb_annotations():
    """EB SIRIUS annotations (2,860 rows) + per-feature m/z & RT."""
    import csv as _csv
    annot = {}
    with open(EB_SIRIUS_ANNOT) as f:
        cols = _csv.DictReader(f, delimiter="\t")
        for row in cols:
            annot[int(row["row ID"])] = row
    mzrt = {}
    with open(EB_TARGETS_CSV) as f:
        for row in _csv.DictReader(f):
            mzrt[int(row["row ID"])] = (float(row["row m/z"]), float(row["row retention time"]))
    missing = [i for i in annot if i not in mzrt]
    if missing:
        sys.exit(f"{len(missing)} annotated EB features missing from {EB_TARGETS_CSV.name}")
    return annot, mzrt


def run_transfer(ppm, rt_tol, cos_min, frag_tol):
    print("parsing spectra...", file=sys.stderr)
    eb_spec = parse_mgf(EB_TARGETS_MGF)
    cur_spec = parse_mgf(FEAT_MGF)
    print(f"  EB spectra: {len(eb_spec)}; this-project spectra with peaks: {len(cur_spec)}",
          file=sys.stderr)

    eb_annot, eb_mzrt = load_eb_annotations()
    fids, fmz, frt, fms2 = load_current_features()
    print(f"  EB SIRIUS-annotated features: {len(eb_annot)}; "
          f"this-project features: {len(fids)} ({fms2.sum()} with MS2)", file=sys.stderr)

    rows = []
    n_assigned = 0
    status_counts: dict[str, int] = {}
    class_counts: dict[str, int] = {}

    for eid, ann in eb_annot.items():
        emz, ert = eb_mzrt[eid]
        cand = np.where(np.abs(fmz - emz) / emz * 1e6 <= ppm)[0]
        cand = cand[np.abs(frt[cand] - ert) <= rt_tol]
        n_cand = int(cand.size)

        scored = []  # (feature_id, cosine)
        for c in cand:
            cid = int(fids[c])
            if cid in cur_spec:
                scored.append((cid, ms2_cosine(eb_spec[eid], cur_spec[cid], frag_tol)))

        # Base row carries the EB chemistry + provenance; resolved below.
        base = {col: ann.get(col, "") for col in CHEM_COLS}
        base["row ID"] = np.nan
        base["sirius_source_feature_id"] = eid
        base["sirius_source_run"] = EB_SOURCE_RUN
        base["annotation_origin"] = "transferred"
        base["n_candidates"] = n_cand
        base["source_mz"] = emz
        base["source_rt"] = ert
        base["feature_mz"] = np.nan
        base["feature_rt"] = np.nan
        base["ppm_error"] = np.nan
        base["rt_delta_min"] = np.nan
        base["ms2_cosine"] = np.nan
        base["match_class"] = ""

        if cand.size == 0:
            base["match_status"] = "no_match"
        elif not scored:
            base["match_status"] = "no_ms2_candidate"
        else:
            scored.sort(key=lambda t: -t[1])
            best_cid, best_cos = scored[0]
            if best_cos < cos_min:
                base["match_status"] = "low_cosine"
            else:
                c_idx = int(np.where(fids == best_cid)[0][0])
                base["match_status"] = "assigned"
                base["row ID"] = best_cid
                base["feature_mz"] = float(fmz[c_idx])
                base["feature_rt"] = float(frt[c_idx])
                base["ppm_error"] = round((float(fmz[c_idx]) - emz) / emz * 1e6, 3)
                base["rt_delta_min"] = round(float(frt[c_idx]) - ert, 3)
                base["ms2_cosine"] = round(best_cos, 4)
                if len(scored) == 1:
                    base["match_class"] = "unique"
                elif scored[1][1] >= cos_min:
                    base["match_class"] = "ms2_tie"
                else:
                    base["match_class"] = "ms2_winner"
                n_assigned += 1

        status_counts[base["match_status"]] = status_counts.get(base["match_status"], 0) + 1
        if base["match_class"]:
            class_counts[base["match_class"]] = class_counts.get(base["match_class"], 0) + 1
        rows.append(base)

    out = pd.DataFrame(rows)
    return out, n_assigned, status_counts, class_counts


# --------------------------------------------------------------------------- #
# Accumulation / merge (native wins over transferred; then structure > confidence)
# --------------------------------------------------------------------------- #
def accumulate(existing: pd.DataFrame | None, new: list[pd.DataFrame]) -> tuple[pd.DataFrame, int, int, int, dict, int]:
    frames = ([existing] if existing is not None and len(existing) else []) + new
    combined = pd.concat(frames, ignore_index=True)

    # Every contributing source is preserved regardless of which row "wins"
    # below -- semicolon-joined, de-duplicated, order-preserving.
    def join_sources(s: pd.Series) -> str:
        seen = []
        for v in s:
            for part in str(v).split(";"):
                if part and part not in seen:
                    seen.append(part)
        return ";".join(seen)

    all_sources = combined.groupby("row ID")["sirius_source_run"].apply(join_sources)
    # all origins contributing to a feature id (transferred / native)
    all_origins = combined.groupby("row ID")["annotation_origin"].apply(join_sources)
    # how many distinct source-run feature ids land on this feature id, which,
    # and the set of formulas those hits carried
    n_sirius_hits = combined.groupby("row ID")["sirius_source_feature_id"].nunique()
    sirius_hit_ids = combined.groupby("row ID")["sirius_source_feature_id"].apply(
        lambda s: ";".join(str(int(x)) for x in sorted(set(int(x) for x in s))))
    sirius_hit_formulas = combined.groupby("row ID")["sirius_formula"].apply(join_sources)
    n_sirius_formulas = combined.groupby("row ID")["sirius_formula"].nunique(dropna=True)
    # conflicts: >1 source hit with DIFFERENT molecular formulas collapsed onto
    # one local feature -- a co-isolated / isobaric pair merged by the
    # Everything-Bagel alignment. The row keeps the highest-priority hit; the
    # rest stay visible in the transfer map. (Concordant duplicate hits, same
    # formula, are tracked via n_sirius_hits but are not conflicts.)
    merged_conflict = n_sirius_formulas > 1

    # precedence: native (0) before transferred (1); then structure; then confidence
    origin_rank = combined["annotation_origin"].map({"native": 0, "transferred": 1}).fillna(2)
    has_structure = (combined["sirius_structure_name"].notna()
                     if "sirius_structure_name" in combined.columns
                     else pd.Series(False, index=combined.index))
    confidence = (combined.get("sirius_structure_confidence", pd.Series(dtype=float))
                  .fillna(-1))
    combined = combined.assign(_origin=origin_rank, _has_structure=has_structure,
                               _confidence=confidence)
    combined = combined.sort_values(["_origin", "_has_structure", "_confidence"],
                                    ascending=[True, False, False])
    winners = (combined.drop_duplicates("row ID", keep="first")
               .drop(columns=["_origin", "_has_structure", "_confidence"]))

    keep = (combined["match_status"].isin(["assigned", "native"]))
    n_assigned = int(keep.sum())
    n_formula = int(combined.loc[keep, "sirius_formula"].notna().sum())
    n_structure = int(combined.loc[keep, "sirius_structure_name"].notna().sum())

    winners = winners.set_index("row ID")
    winners["sirius_source_run"] = all_sources
    winners["annotation_origin"] = all_origins
    winners["n_sirius_hits"] = n_sirius_hits
    winners["sirius_hit_ids"] = sirius_hit_ids
    winners["sirius_hit_formulas"] = sirius_hit_formulas
    winners["merged_conflict"] = merged_conflict
    winners = winners.reset_index().sort_values("row ID")

    wins = winners[winners["match_status"].isin(["assigned", "native"])]
    per_source = wins.groupby("annotation_origin").size().to_dict()
    n_conflict = int(winners["merged_conflict"].sum())
    return winners, n_assigned, n_formula, n_structure, per_source, n_conflict


# --------------------------------------------------------------------------- #
# main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ppm", type=float, default=10.0,
                    help="precursor m/z tolerance (ppm) for candidate search [10]")
    ap.add_argument("--rt", type=float, default=0.5,
                    help="retention-time tolerance (min) for candidate search [0.5]")
    ap.add_argument("--cos-min", type=float, default=0.7,
                    help="minimum MS2 cosine to accept an assignment [0.7]")
    ap.add_argument("--frag-tol", type=float, default=0.05,
                    help="fragment m/z tolerance (Da) for the cosine [0.05]")
    ap.add_argument("--native-merged", type=Path, default=None,
                    help="merged dir of a COMPLETED native SIRIUS run on this "
                         "project's aligned_features_filled.mgf, to merge in "
                         "(native rows win over transferred for the same feature)")
    ap.add_argument("--native-label", default=None,
                    help="source_run label for the native run (default: "
                         "--native-merged absolute path)")
    ap.add_argument("--fresh", action="store_true",
                    help="ignore existing sirius_annotations.tsv and rebuild")
    args = ap.parse_args()

    # existing accumulated table (from prior transfer and/or prior native merge)
    existing = None
    if not args.fresh and OUT_ANNOT.exists():
        existing = pd.read_csv(OUT_ANNOT, sep="\t")
        # guard against loading a file written before provenance columns existed
        for col in PROV_COLS + CHEM_COLS:
            if col not in existing.columns:
                existing[col] = np.nan if col not in CHEM_COLS else ""
        print(f"accumulating onto {len(existing)} previously-imported features "
              f"({OUT_ANNOT.name})", file=sys.stderr)

    new_frames = []

    # 1) the transfer itself (unless only patch-merging a native run)
    if not args.native_merged or not (existing is not None and not existing.empty):
        transfer_df, n_assigned, status_counts, class_counts = run_transfer(
            args.ppm, args.rt, args.cos_min, args.frag_tol)
        print("transfer assignment counts:", dict(sorted(status_counts.items())),
              file=sys.stderr)
        print(f"  assigned {n_assigned}/{len(transfer_df)} "
              f"({n_assigned / len(transfer_df) * 100:.1f}%) "
              f"classes={dict(sorted(class_counts.items()))}", file=sys.stderr)
        OUT_MAP.parent.mkdir(parents=True, exist_ok=True)
        transfer_df.sort_values("sirius_source_feature_id").to_csv(OUT_MAP, sep="\t", index=False)
        # only assigned rows carry an annotation worth accumulating
        transfer_assigned = transfer_df[transfer_df["match_status"] == "assigned"].copy()
        print(f"  wrote {OUT_MAP.name} ({len(transfer_df)} EB features, "
              f"{len(transfer_assigned)} assigned)", file=sys.stderr)
        new_frames.append(transfer_assigned)
    else:
        print("skipping transfer (existing table present and --native-merged given)",
              file=sys.stderr)

    # 2) optionally merge a native run
    if args.native_merged is not None:
        label = args.native_label or str(args.native_merged.resolve())
        native_df = distill_native_run(args.native_merged, label)
        print(f"[{label}] native: {len(native_df)} features, "
              f"{native_df['sirius_formula'].notna().sum()} formula, "
              f"{native_df['sirius_structure_name'].notna().sum()} structure",
              file=sys.stderr)
        for col in PROV_COLS + CHEM_COLS:
            if col not in native_df.columns:
                native_df[col] = np.nan if col not in CHEM_COLS else ""
        new_frames.append(native_df)

    if not new_frames:
        sys.exit("nothing to do: pass a transfer config or --native-merged")

    winners, n_assigned, n_formula, n_structure, per_source, n_conflict = accumulate(existing, new_frames)
    # drop never-assigned rows from the accumulated annotation table
    winners = winners[winners["match_status"].isin(["assigned", "native"])]

    OUT_ANNOT.parent.mkdir(parents=True, exist_ok=True)
    winners.sort_values("row ID").to_csv(OUT_ANNOT, sep="\t", index=False)
    print(f"accumulated {len(winners)} assigned feature ids -> {OUT_ANNOT}",
          file=sys.stderr)
    print(f"  formula {n_formula}; structure {n_structure}; by origin {per_source}",
          file=sys.stderr)
    print(f"  merged-conflict features (n_sirius_hits>1): {n_conflict}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
