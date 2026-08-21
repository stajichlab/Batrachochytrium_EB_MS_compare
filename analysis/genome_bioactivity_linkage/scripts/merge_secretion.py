"""Combine SignalP, DeepTMHMM, and PredGPI calls into a predicted-extracellular
protein set (spec Stage 2)."""
import gzip
from pathlib import Path

import pandas as pd

from parse_deeptmhmm import has_tm_helix_outside_signal


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_signalp_gff3(path: Path) -> pd.DataFrame:
    """Load a real BFD SignalP GFF3 file into ``protein_id,
    is_signal_peptide, cleavage_site`` (one row per protein).

    Real SignalP GFF3 (verified against
    ``.../results/function/signalp/9f/FF9F6419.signalp.gff3.gz``) starts
    with a ``## gff-version 3`` line, then one feature line per predicted
    signal peptide of type ``signal_peptide`` -- only SignalP-POSITIVE
    proteins appear in the file at all (there are no explicit negative
    rows). Column 1 (seqid) is ``"<protein_id> <protein_id>"`` (the
    protein id repeated, space-separated) -- split on whitespace and take
    the first token. The cleavage site is column 5 (the feature's end
    coordinate).
    """
    rows = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 5 or fields[2] != "signal_peptide":
                continue
            protein_id = fields[0].split()[0]
            cleavage_site = int(fields[4])
            rows.append(
                {
                    "protein_id": protein_id,
                    "is_signal_peptide": True,
                    "cleavage_site": cleavage_site,
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "is_signal_peptide", "cleavage_site"])


def load_predgpi_gff3(path: Path) -> pd.DataFrame:
    """Load a real BFD PredGPI GFF3 file into ``protein_id,
    has_gpi_anchor`` (one row per protein mentioned anywhere in the file).

    Real PredGPI GFF3 (verified against
    ``.../results/function/predgpi/9c/FE87067D.predgpi.gff3.gz``) has no
    header/comment lines; the feature-type column (column 3) is either
    ``Chain`` (no GPI anchor -- the vast majority of rows) or
    ``GPI-anchor`` (has a GPI anchor). There is no explicit boolean column
    -- ``has_gpi_anchor`` is derived from whether a protein has any
    ``GPI-anchor``-type feature row.
    """
    gpi_anchor_proteins = set()
    all_proteins = set()
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            protein_id = fields[0].split()[0]
            all_proteins.add(protein_id)
            if fields[2] == "GPI-anchor":
                gpi_anchor_proteins.add(protein_id)
    rows = [
        {"protein_id": protein_id, "has_gpi_anchor": protein_id in gpi_anchor_proteins}
        for protein_id in sorted(all_proteins)
    ]
    return pd.DataFrame(rows, columns=["protein_id", "has_gpi_anchor"])


def predicted_extracellular(
    signalp: pd.DataFrame, deeptmhmm_gff3: pd.DataFrame, predgpi: pd.DataFrame
) -> pd.DataFrame:
    merged = signalp.merge(predgpi, on="protein_id", how="left")
    merged["has_gpi_anchor"] = merged["has_gpi_anchor"].fillna(False).astype(bool)

    def _disqualifying_tm(row):
        if not row["is_signal_peptide"]:
            return False  # irrelevant once excluded by signalp_positive below
        return has_tm_helix_outside_signal(deeptmhmm_gff3, row["protein_id"], row["cleavage_site"])

    merged["has_disqualifying_tm"] = merged.apply(_disqualifying_tm, axis=1).astype(bool)
    merged["signalp_positive"] = merged["is_signal_peptide"]
    merged["signal_cleavage_site"] = merged["cleavage_site"]
    merged["is_extracellular"] = (
        merged["signalp_positive"]
        & ~merged["has_disqualifying_tm"]
        & ~merged["has_gpi_anchor"]
    )
    return merged[
        [
            "protein_id",
            "signalp_positive",
            "signal_cleavage_site",
            "has_disqualifying_tm",
            "has_gpi_anchor",
            "is_extracellular",
        ]
    ]
