import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_linkage_tables  # noqa: E402
from build_linkage_tables import build_gene_domain_table, compute_bgc_proteins  # noqa: E402


def test_build_gene_domain_table_merges_pfam_bgc_and_secretion_and_crossref():
    pfam_calls = pd.DataFrame(
        [{"protein_id": "protA", "pfam_id": "PF01397", "domain_name": "TPS1", "evalue": 1e-40, "score": 140.0, "family": "terpene_synthase"}]
    )
    bgc_proteins = {"protA"}
    secretion = pd.DataFrame(
        [{"protein_id": "protA", "signalp_positive": True, "signal_cleavage_site": 20,
          "has_disqualifying_tm": False, "has_gpi_anchor": False, "is_extracellular": True}]
    )
    rbh_confirmed_proteins = {"protA"}

    result = build_gene_domain_table(pfam_calls, bgc_proteins, secretion, rbh_confirmed_proteins)
    row = result[result["protein_id"] == "protA"].iloc[0]
    assert row["family"] == "terpene_synthase"
    assert bool(row["has_bgc_context"]) is True
    assert bool(row["is_cross_ref_confirmed"]) is True
    assert bool(row["is_extracellular"]) is True


def test_build_gene_domain_table_defaults_false_when_missing_from_bgc_or_rbh_sets():
    pfam_calls = pd.DataFrame(
        [{"protein_id": "protZ", "pfam_id": "PF00109", "domain_name": "ketoacyl-synt", "evalue": 1e-30, "score": 100.0, "family": "pks"}]
    )
    secretion = pd.DataFrame(
        [{"protein_id": "protZ", "signalp_positive": True, "signal_cleavage_site": 20,
          "has_disqualifying_tm": False, "has_gpi_anchor": False, "is_extracellular": True}]
    )
    result = build_gene_domain_table(pfam_calls, bgc_proteins=set(), secretion=secretion, rbh_confirmed_proteins=set())
    row = result[result["protein_id"] == "protZ"].iloc[0]
    assert bool(row["has_bgc_context"]) is False
    assert bool(row["is_cross_ref_confirmed"]) is False


def test_compute_bgc_proteins_uses_real_genomic_coords_not_fullhmmer_membership():
    # protA's genomic coords overlap the BGC region -> included.
    # protB has a PFAM call but its genomic coords fall OUTSIDE any BGC
    # region -- it must NOT be flagged as in-BGC even though a naive
    # "protein appears in fullhmmer hits" shortcut would wrongly include
    # nearly every domain hit (see module docstring deviation 1).
    pfam_calls = pd.DataFrame(
        [
            {"protein_id": "protA", "pfam_id": "PF01397", "family": "terpene_synthase"},
            {"protein_id": "protB", "pfam_id": "PF00109", "family": "pks"},
            {"protein_id": "protC", "pfam_id": "PF00067", "family": "p450"},  # no coords available
        ]
    )
    protein_coords = {
        "protA": ("contig1", 1200, 2200),
        "protB": ("contig1", 6000, 6300),
    }
    regions = [{"record_id": "contig1", "start": 1000, "end": 5000}]

    bgc_proteins = compute_bgc_proteins(pfam_calls, protein_coords, regions)

    assert bgc_proteins == {"protA"}


def test_load_compound_log2fc_qvalue_keeps_lowest_qvalue_liq_row_and_excludes_others(tmp_path, monkeypatch):
    fixture = tmp_path / "all_significant_features_summary.tsv"
    fixture.write_text(
        "row_id\tspecies\tcomparison\tlog2FC_a_over_b\tq_value\n"
        # row 1: two liq contrasts for dendrobatidis -- keep the lower q_value one
        "1\tdendrobatidis\tdendrobatidis_liq_Developed_vs_spore_Developed\t2.0\t0.05\n"
        "1\tdendrobatidis\tdendrobatidis_liq_Zoospore_vs_liq_Developed\t-3.0\t0.01\n"
        # row 2: a non-liq contrast only -- must be excluded (no liq evidence)
        "2\tdendrobatidis\tdendrobatidis_spore_Zoospore_vs_spore_Developed\t4.0\t0.001\n"
        # row 3: different species -- must be excluded when querying dendrobatidis
        "3\tsalamandrivorans\tsalamandrivorans_liq_Developed_vs_spore_Developed\t1.0\t0.02\n"
    )
    monkeypatch.setattr(build_linkage_tables, "DIFFERENTIAL_FEATURES_PATH", fixture)

    result = build_linkage_tables.load_compound_log2fc_qvalue("dendrobatidis")

    assert set(result["row_id"]) == {1}
    row = result[result["row_id"] == 1].iloc[0]
    assert row["log2fc"] == -3.0
    assert row["q_value"] == 0.01
