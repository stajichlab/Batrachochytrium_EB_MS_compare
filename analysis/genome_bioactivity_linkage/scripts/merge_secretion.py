"""Combine SignalP, DeepTMHMM, and PredGPI calls into a predicted-extracellular
protein set (spec Stage 2)."""
import pandas as pd

from parse_deeptmhmm import has_tm_helix_outside_signal


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
