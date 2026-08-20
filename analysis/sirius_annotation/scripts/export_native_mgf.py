#!/usr/bin/env python3
"""
Export native SIRIUS targets MGF from the Everything-Bagel feature MGF.
Pulls the MS2 block for each target row ID out of
data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/aligned_features_filled.mgf
into analysis/sirius_annotation/sirius_native_targets.mgf.

Join key: SCANS (== FEATURE_ID, == aligned feature "row ID") in the feature MGF
equals "row ID" in sirius_native_targets.csv -- verified: block headers carry
TITLE=SCAN=N / FEATURE_ID=N / SCANS=N with the aligned feature id.

Validation (blocks that fail are DROPPED and reported -- SIRIUS cannot use them):
  * CHARGE tag present and == the expected target charge (default 1+)
  * PEPMASS > 0 (the feature MGF has degenerate blocks with PEPMASS=0.0)
  * at least one peak with intensity > 0

Usage:
    python3 scripts/export_native_mgf.py [--charge 1]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent.parent.parent
TARGETS_CSV = REPO / "analysis" / "sirius_annotation" / "sirius_native_targets.csv"
SOURCE_MGF = (
    REPO
    / "data"
    / "raw"
    / "gnps2_e9838293_bagel"
    / "nf_output"
    / "feature_finding"
    / "aligned_features_filled.mgf"
)
OUT_MGF = REPO / "analysis" / "sirius_annotation" / "sirius_native_targets.mgf"


def parse_mgf_blocks(text: str):
    """Yield (scans_id, full_block_text) for every BEGIN IONS...END IONS block."""
    block_lines = []
    scans = None
    for line in text.splitlines(keepends=True):
        if line.startswith("BEGIN IONS"):
            block_lines = [line]
            scans = None
        elif line.startswith("END IONS"):
            block_lines.append(line)
            yield scans, "".join(block_lines)
        else:
            block_lines.append(line)
            if line.startswith("SCANS="):
                scans = int(line.strip().split("=", 1)[1])


def _pepmass_of(block: str) -> float:
    for line in block.splitlines():
        if line.startswith("PEPMASS="):
            try:
                return float(line.split("=", 1)[1].split()[0])
            except (ValueError, IndexError):
                return 0.0
    return 0.0


def _charge_of(block: str) -> str | None:
    for line in block.splitlines():
        if line.startswith("CHARGE="):
            return line.split("=", 1)[1].strip()
    return None


def _has_positive_peak(block: str) -> bool:
    for line in block.splitlines():
        parts = line.replace(",", " ").split()
        if len(parts) == 2:
            try:
                if float(parts[1]) > 0.0:
                    return True
            except ValueError:
                continue
    return False


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--charge", type=int, default=1, help="Expected target charge (default 1)")
    args = ap.parse_args()
    expected_charge = f"{args.charge}+"

    targets = pd.read_csv(TARGETS_CSV)
    target_ids = sorted(set(targets["row ID"].astype(int)))
    print(f"looking up {len(target_ids)} target row IDs in {SOURCE_MGF.name}", file=sys.stderr)

    with SOURCE_MGF.open() as fh:
        text = fh.read()

    found = {}
    bad_charge = {}
    for scans, block in parse_mgf_blocks(text):
        if scans not in target_ids:
            continue
        if _charge_of(block) != expected_charge:
            bad_charge[scans] = f"charge={_charge_of(block)}"
            continue
        found[scans] = block

    missing = [t for t in target_ids if t not in found and t not in bad_charge]
    no_pepmass = {s: b for s, b in found.items() if _pepmass_of(b) <= 0.0}
    for s in no_pepmass:
        found.pop(s)
    no_peaks = {s: b for s, b in found.items() if not _has_positive_peak(b)}
    for s in no_peaks:
        found.pop(s)

    dropped = 0
    for label, d in (("missing block", missing), ("bad CHARGE", bad_charge)):
        if d:
            print(f"DROPPED {len(d)} targets ({label}):", file=sys.stderr)
            for s in sorted(d)[:20]:
                print(f"  row ID {s} {d[s] if isinstance(d[s], str) else ''}", file=sys.stderr)
            dropped += len(d)
    for label, d in (("PEPMASS<=0", no_pepmass), ("no positive peaks", no_peaks)):
        if d:
            print(f"DROPPED {len(d)} targets ({label}): {sorted(d)[:20]}", file=sys.stderr)
            dropped += len(d)

    OUT_MGF.parent.mkdir(parents=True, exist_ok=True)
    with OUT_MGF.open("w") as fh:
        for scans in sorted(found):
            fh.write(found[scans])
    print(
        f"wrote {len(found)}/{len(target_ids)} usable spectra (dropped {dropped}) -> {OUT_MGF}",
        file=sys.stderr,
    )
    if len(found) != len(target_ids):
        print(
            f"WARNING: {len(target_ids) - len(found)} target(s) have no usable MS2; "
            "they cannot be SIRIUS-annotated",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
