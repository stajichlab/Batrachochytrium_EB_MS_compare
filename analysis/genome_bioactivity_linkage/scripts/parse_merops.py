"""Parse a MEROPS blastp result (spec: peptidase-family/clan classification,
companion to parse_pfam_domains.py for the biosynthetic-domain side)."""
import gzip
from pathlib import Path

import pandas as pd

MEROPS_FAMILIES_TAB = Path("/srv/projects/db/MEROPS/124/merops_lib.families.tab")

# MEROPS catalytic-type prefix on the family code (e.g. "S08" -> serine).
# "I" is a peptidase-INHIBITOR family, not a peptidase itself -- kept
# distinct so callers can exclude it rather than mis-report it as a hit.
CATALYTIC_TYPE = {
    "S": "serine",
    "C": "cysteine",
    "A": "aspartic",
    "M": "metallo",
    "T": "threonine",
    "G": "glutamic",
    "N": "asparagine",
    "U": "unknown",
    "I": "inhibitor",
}


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_merops_blasttab(path: Path) -> pd.DataFrame:
    """Load a real MEROPS blastp ``-outfmt 6`` result (protein_id, MER
    subject id, pident, length, mismatch, gapopen, qstart, qend, sstart,
    send, evalue, bitscore) -- standard BLAST tabular columns, no header.
    """
    cols = [
        "protein_id", "mer_id", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore",
    ]
    rows = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if len(fields) < 12:
                continue
            rows.append(fields[:12])
    df = pd.DataFrame(rows, columns=cols)
    if df.empty:
        return df
    for c in ("pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore"):
        df[c] = pd.to_numeric(df[c])
    return df


def load_merops_families(path: Path = MEROPS_FAMILIES_TAB) -> pd.DataFrame:
    """Load the MEROPS id -> clan/family lookup table (``MER_id\\tclan\\tfamily``,
    e.g. ``MER0000002\\tS01A\\tS01.001``)."""
    df = pd.read_csv(path, sep="\t", names=["mer_id", "clan", "family"])
    df["catalytic_type"] = df["family"].str[0].map(CATALYTIC_TYPE).fillna("unknown")
    return df


def best_merops_hit(blasttab: pd.DataFrame, families: pd.DataFrame) -> pd.DataFrame:
    """Reduce a blasttab to one best (lowest e-value) MEROPS hit per protein_id,
    annotated with family/clan/catalytic_type."""
    if blasttab.empty:
        return pd.DataFrame(
            columns=["protein_id", "mer_id", "pident", "evalue", "bitscore", "clan", "family", "catalytic_type"]
        )
    best = blasttab.sort_values("evalue", ascending=True).drop_duplicates("protein_id")
    merged = best.merge(families, on="mer_id", how="left")
    return merged[["protein_id", "mer_id", "pident", "evalue", "bitscore", "clan", "family", "catalytic_type"]]
