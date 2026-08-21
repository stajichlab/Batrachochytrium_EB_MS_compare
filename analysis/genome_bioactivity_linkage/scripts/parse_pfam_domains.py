"""Parse BFD's hmmsearch --cut_ga domtblout output and classify Pfam hits into
biosynthetic domain families (see domain_families.py)."""
import gzip
from pathlib import Path

import pandas as pd

from domain_families import classify_pfam


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_domtblout(path: Path) -> pd.DataFrame:
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(None, 22)
            protein_id = fields[0]
            domain_name = fields[3]
            pfam_id_versioned = fields[4]
            evalue_full = float(fields[6])
            score_full = float(fields[7])
            rows.append(
                {
                    "protein_id": protein_id,
                    "pfam_id": pfam_id_versioned.split(".")[0],
                    "domain_name": domain_name,
                    "evalue": evalue_full,
                    "score": score_full,
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "pfam_id", "domain_name", "evalue", "score"])


def classify_domains(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["family"] = out["pfam_id"].map(classify_pfam)
    return out.dropna(subset=["family"]).reset_index(drop=True)
