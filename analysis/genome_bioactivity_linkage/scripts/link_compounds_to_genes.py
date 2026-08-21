"""Final tiered candidate-gene table for each background-filtered, class-mapped
liquid-fraction compound (spec Stage 4). Tiered/lexicographic ranking only --
never a single weighted composite score."""
import pandas as pd

from domain_families import COMPOUND_CLASS_TO_FAMILY


def assign_tier(has_bgc_context: bool, is_cross_ref_confirmed: bool) -> int:
    if has_bgc_context and is_cross_ref_confirmed:
        return 1
    if has_bgc_context or is_cross_ref_confirmed:
        return 2
    return 3


def build_candidate_table(compounds: pd.DataFrame, gene_domains: pd.DataFrame) -> pd.DataFrame:
    extracellular = gene_domains[gene_domains["is_extracellular"]]
    rows = []
    for _, compound in compounds.iterrows():
        family = COMPOUND_CLASS_TO_FAMILY.get(compound["compound_class"])
        if family is None:
            continue  # unmapped compound class -- left out per spec, not force-fit
        matches = extracellular[extracellular["family"] == family]
        for _, gene in matches.iterrows():
            tier = assign_tier(gene["has_bgc_context"], gene["is_cross_ref_confirmed"])
            rows.append(
                {
                    "compound_row_id": compound["row_id"],
                    "compound_class": compound["compound_class"],
                    "candidate_protein_id": gene["protein_id"],
                    "domain_family": family,
                    "tier": tier,
                    "has_bgc_context": gene["has_bgc_context"],
                    "is_cross_ref_confirmed": gene["is_cross_ref_confirmed"],
                    "is_extracellular": gene["is_extracellular"],
                    "compound_log2fc": compound["log2fc"],
                    "compound_q_value": compound["q_value"],
                }
            )
    table = pd.DataFrame(
        rows,
        columns=[
            "compound_row_id", "compound_class", "candidate_protein_id", "domain_family",
            "tier", "has_bgc_context", "is_cross_ref_confirmed", "is_extracellular",
            "compound_log2fc", "compound_q_value",
        ],
    )
    if table.empty:
        return table
    # Ranked tier-first, not compound-first: ``compound_log2fc`` is constant
    # within every compound_row_id group (it's a compound-level attribute,
    # not a per-row value), so a compound-first sort would make the
    # fold-change key inert as a tie-breaker and the table wouldn't be
    # globally ranked by tier the way the spec describes ("tie-breaker
    # within a tier" implies tier-first ranking). compound_row_id is kept as
    # the final key purely for deterministic, stable ordering among
    # otherwise-tied (tier, |fc|) rows.
    table["_abs_fc"] = table["compound_log2fc"].abs()
    table = table.sort_values(
        ["tier", "_abs_fc", "compound_row_id"], ascending=[True, False, True]
    ).drop(columns="_abs_fc").reset_index(drop=True)
    return table
