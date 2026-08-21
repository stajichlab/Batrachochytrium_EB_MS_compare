import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from link_compounds_to_genes import assign_tier, build_candidate_table  # noqa: E402


def test_assign_tier():
    assert assign_tier(has_bgc_context=True, is_cross_ref_confirmed=True) == 1
    assert assign_tier(has_bgc_context=False, is_cross_ref_confirmed=True) == 2
    assert assign_tier(has_bgc_context=True, is_cross_ref_confirmed=False) == 2
    assert assign_tier(has_bgc_context=False, is_cross_ref_confirmed=False) == 3


def test_build_candidate_table_only_matches_extracellular_matching_family():
    compounds = pd.DataFrame(
        [
            {"row_id": 100, "compound_class": "terpenoid", "log2fc": 3.0, "q_value": 0.01},
            {"row_id": 200, "compound_class": "polyketide", "log2fc": 1.0, "q_value": 0.2},
        ]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "protA",
                "family": "terpene_synthase",
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
            {
                "protein_id": "protB",
                "family": "terpene_synthase",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": False,  # excluded: not secreted
            },
            {
                "protein_id": "protC",
                "family": "pks",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
        ]
    )
    table = build_candidate_table(compounds, gene_domains)
    assert set(table["candidate_protein_id"]) == {"protA", "protC"}
    row_a = table[table["candidate_protein_id"] == "protA"].iloc[0]
    assert row_a["compound_row_id"] == 100
    assert row_a["tier"] == 1
    row_c = table[table["candidate_protein_id"] == "protC"].iloc[0]
    assert row_c["compound_row_id"] == 200
    assert row_c["tier"] == 2


def test_build_candidate_table_sorted_by_compound_then_tier_then_fc():
    compounds = pd.DataFrame(
        [{"row_id": 1, "compound_class": "terpenoid", "log2fc": 5.0, "q_value": 0.01}]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "low_tier",
                "family": "terpene_synthase",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": True,
            },
            {
                "protein_id": "high_tier",
                "family": "terpene_synthase",
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
        ]
    )
    table = build_candidate_table(compounds, gene_domains)
    assert list(table["candidate_protein_id"]) == ["high_tier", "low_tier"]
