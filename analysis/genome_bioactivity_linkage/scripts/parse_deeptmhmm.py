"""Parse DeepTMHMM TMRs.gff3 output and apply the SignalP/DeepTMHMM overlap rule.

Rule (spec Stage 2): a TM helix call that falls entirely within the
SignalP-cleaved N-terminal region does NOT disqualify a protein from being
"predicted extracellular"; only a TM helix that extends into the mature
chain (past the cleavage site) does.
"""
from pathlib import Path

import pandas as pd


def parse_tmrs_gff3(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            protein_id, region_type, start, end = line.split("\t")
            rows.append(
                {
                    "protein_id": protein_id,
                    "region_type": region_type,
                    "start": int(start),
                    "end": int(end),
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "region_type", "start", "end"])


def has_tm_helix_outside_signal(
    df: pd.DataFrame, protein_id: str, signal_cleavage_site: int | None
) -> bool:
    tm_rows = df[(df["protein_id"] == protein_id) & (df["region_type"] == "TMhelix")]
    if tm_rows.empty:
        return False
    if signal_cleavage_site is None:
        return True
    return bool((tm_rows["end"] > signal_cleavage_site).any())
