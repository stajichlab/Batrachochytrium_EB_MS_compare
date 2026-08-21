"""Reduce two one-directional DIAMOND blastp -outfmt 6 hit lists to reciprocal best hits."""
from pathlib import Path

import pandas as pd

_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def _top_hits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=_COLUMNS)
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid", keep="first")
    return df.set_index("qseqid")


def reciprocal_best_hits(fwd_tsv: Path, rev_tsv: Path) -> pd.DataFrame:
    fwd_top = _top_hits(fwd_tsv)  # bfd_id -> best reference hit
    rev_top = _top_hits(rev_tsv)  # reference_id -> best bfd hit

    rows = []
    for bfd_id, fwd_row in fwd_top.iterrows():
        ref_id = fwd_row["sseqid"]
        if ref_id in rev_top.index and rev_top.loc[ref_id, "sseqid"] == bfd_id:
            rows.append(
                {
                    "bfd_protein_id": bfd_id,
                    "reference_protein_id": ref_id,
                    "pident": fwd_row["pident"],
                    "evalue": fwd_row["evalue"],
                }
            )
    return pd.DataFrame(rows, columns=["bfd_protein_id", "reference_protein_id", "pident", "evalue"])
