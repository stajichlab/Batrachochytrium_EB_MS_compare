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
            {"row_id": 100, "compound_class": "Terpenoids", "log2fc": 3.0, "q_value": 0.01},
            {"row_id": 200, "compound_class": "Polyketides", "log2fc": 1.0, "q_value": 0.2},
        ]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "protA",
                "family": "terpene_synthase",
                "n_domain_hits": 1,
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
            {
                "protein_id": "protB",
                "family": "terpene_synthase",
                "n_domain_hits": 1,
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": False,  # excluded: not secreted
            },
            {
                "protein_id": "protC",
                "family": "pks",
                "n_domain_hits": 2,
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
    assert row_a["n_domain_hits"] == 1
    row_c = table[table["candidate_protein_id"] == "protC"].iloc[0]
    assert row_c["compound_row_id"] == 200
    assert row_c["tier"] == 2
    # n_domain_hits (evidence preserved from build_gene_domain_table's I1
    # dedup) must actually reach the final output table, not just be
    # computed and dropped -- see GENOME_BIOACTIVITY_LINKAGE.md residual #3.
    assert row_c["n_domain_hits"] == 2


def test_build_candidate_table_sorted_by_tier_then_fc_then_compound():
    # Two different compounds (genuinely different fold-changes) both
    # matching "terpene_synthase", against genes spanning all three tiers,
    # to prove the table is ranked tier-first (not compound-first) with
    # |compound_log2fc| breaking ties within a tier -- see
    # GENOME_BIOACTIVITY_LINKAGE.md I2.
    compounds = pd.DataFrame(
        [
            {"row_id": 2, "compound_class": "Terpenoids", "log2fc": 1.0, "q_value": 0.03},
            {"row_id": 1, "compound_class": "Terpenoids", "log2fc": 5.0, "q_value": 0.01},
        ]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "low_tier",
                "family": "terpene_synthase",
                "n_domain_hits": 1,
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": True,
            },
            {
                "protein_id": "high_tier",
                "family": "terpene_synthase",
                "n_domain_hits": 1,
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
        ]
    )
    table = build_candidate_table(compounds, gene_domains)
    # Tier-first: both "high_tier" rows (tier 1) precede both "low_tier"
    # rows (tier 3), regardless of which compound they're paired with.
    assert list(table["tier"]) == [1, 1, 3, 3]
    assert list(table["candidate_protein_id"]) == ["high_tier", "high_tier", "low_tier", "low_tier"]
    # Within tier 1 (a tie on protein), the higher |fc| compound (row_id=1,
    # fc=5.0) ranks before the lower |fc| compound (row_id=2, fc=1.0) --
    # proving the fold-change tie-breaker is not inert.
    tier1 = table[table["tier"] == 1]
    assert list(tier1["compound_row_id"]) == [1, 2]
    tier3 = table[table["tier"] == 3]
    assert list(tier3["compound_row_id"]) == [1, 2]
