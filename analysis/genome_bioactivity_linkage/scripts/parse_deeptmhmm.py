"""Parse DeepTMHMM TMRs.gff3 output and apply the SignalP/DeepTMHMM overlap rule.

Rule (spec Stage 2): a TM helix call that falls entirely within the
SignalP-cleaved N-terminal region does NOT disqualify a protein from being
"predicted extracellular"; only a TM helix that extends into the mature
chain (past the cleavage site) does.
"""
from pathlib import Path

import pandas as pd


def parse_tmrs_gff3(path: Path) -> pd.DataFrame:
    """Parse a real DeepTMHMM ``TMRs.gff3`` file.

    Real DeepTMHMM output (verified against a live run's output, e.g.
    ``/rhome/jstajich/projects/nf/nf_funannotate1/tests/output/
    deeptmhmm_gpu_test/27596677/TMRs.gff3``) has two quirks a naive GFF3
    parser doesn't expect:

    1. Data lines are tab-padded to 8 fields (e.g.
       ``P40231.2\\tinside\\t1\\t332\\t\\t\\t\\t``) -- only the first 4 fields
       (protein_id, region_type, start, end) are meaningful.
    2. Records are separated by lines containing just ``//`` (not a ``#``
       comment), in addition to the ``# ...`` header/comment lines.
    """
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#") or line == "//":
                continue
            protein_id, region_type, start, end = line.split("\t")[:4]
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
