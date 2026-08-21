import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import build_linkage_tables  # noqa: E402
from build_linkage_tables import (  # noqa: E402
    _strip_transcript_suffix,
    build_gene_domain_table,
    compute_bgc_proteins,
)


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
    assert row["n_domain_hits"] == 1
    assert bool(row["has_bgc_context"]) is True
    assert bool(row["is_cross_ref_confirmed"]) is True
    assert bool(row["is_extracellular"]) is True


def test_build_gene_domain_table_dedups_multiple_domain_hits_same_family():
    # A protein with 2 domain hits in the same family (e.g. a PKS with both
    # a KS domain PF00109 and an AT domain PF00698, both -> "pks") must
    # collapse to exactly ONE output row, with n_domain_hits == 2 rather
    # than losing that evidence entirely -- see GENOME_BIOACTIVITY_LINKAGE.md I1.
    pfam_calls = pd.DataFrame(
        [
            {"protein_id": "protPKS", "pfam_id": "PF00109", "domain_name": "ketoacyl-synt", "evalue": 1e-30, "score": 100.0, "family": "pks"},
            {"protein_id": "protPKS", "pfam_id": "PF00698", "domain_name": "Acyl_transf_1", "evalue": 1e-20, "score": 80.0, "family": "pks"},
        ]
    )
    secretion = pd.DataFrame(
        [{"protein_id": "protPKS", "signalp_positive": True, "signal_cleavage_site": 20,
          "has_disqualifying_tm": False, "has_gpi_anchor": False, "is_extracellular": True}]
    )
    result = build_gene_domain_table(pfam_calls, bgc_proteins=set(), secretion=secretion, rbh_confirmed_proteins=set())
    assert len(result) == 1
    row = result.iloc[0]
    assert row["protein_id"] == "protPKS"
    assert row["family"] == "pks"
    assert row["n_domain_hits"] == 2


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
        # row 1: two genuine liq-vs-spore contrasts for dendrobatidis, both
        # positive (liq-enriched) -- keep the lower q_value one.
        "1\tdendrobatidis\tdendrobatidis_liq_Developed_vs_spore_Developed\t2.0\t0.05\n"
        "1\tdendrobatidis\tdendrobatidis_liq_Zoospore_vs_spore_Zoospore\t3.0\t0.01\n"
        # row 2: a non-liq contrast only -- must be excluded (no liq evidence)
        "2\tdendrobatidis\tdendrobatidis_spore_Zoospore_vs_spore_Developed\t4.0\t0.001\n"
        # row 3: different species -- must be excluded when querying dendrobatidis
        "3\tsalamandrivorans\tsalamandrivorans_liq_Developed_vs_spore_Developed\t1.0\t0.02\n"
        # row 4: liq-vs-spore but NEGATIVE fc (spore-enriched, wrong
        # direction) -- must be excluded even though comparison matches.
        "4\tdendrobatidis\tdendrobatidis_liq_Developed_vs_spore_Developed\t-2.5\t0.001\n"
        # row 5: liq-vs-liq life-stage contrast (contains 'liq' twice, no
        # 'vs_spore') -- must be excluded regardless of fc sign, since it
        # carries no liq-vs-spore enrichment information at all.
        "5\tdendrobatidis\tdendrobatidis_liq_Zoospore_vs_liq_Developed\t5.0\t0.0001\n"
    )
    monkeypatch.setattr(build_linkage_tables, "DIFFERENTIAL_FEATURES_PATH", fixture)

    result = build_linkage_tables.load_compound_log2fc_qvalue("dendrobatidis")

    assert set(result["row_id"]) == {1}
    row = result[result["row_id"] == 1].iloc[0]
    assert row["log2fc"] == 3.0
    assert row["q_value"] == 0.01


def test_load_compound_log2fc_qvalue_liq_vs_spore_positive_fc_kept(tmp_path, monkeypatch):
    fixture = tmp_path / "all_significant_features_summary.tsv"
    fixture.write_text(
        "row_id\tspecies\tcomparison\tlog2FC_a_over_b\tq_value\n"
        "10\tdendrobatidis\tdendrobatidis_liq_Developed_vs_spore_Developed\t1.5\t0.02\n"
    )
    monkeypatch.setattr(build_linkage_tables, "DIFFERENTIAL_FEATURES_PATH", fixture)
    result = build_linkage_tables.load_compound_log2fc_qvalue("dendrobatidis")
    assert set(result["row_id"]) == {10}


def test_load_compound_log2fc_qvalue_liq_vs_spore_negative_fc_excluded(tmp_path, monkeypatch):
    fixture = tmp_path / "all_significant_features_summary.tsv"
    fixture.write_text(
        "row_id\tspecies\tcomparison\tlog2FC_a_over_b\tq_value\n"
        "11\tdendrobatidis\tdendrobatidis_liq_Developed_vs_spore_Developed\t-1.5\t0.02\n"
    )
    monkeypatch.setattr(build_linkage_tables, "DIFFERENTIAL_FEATURES_PATH", fixture)
    result = build_linkage_tables.load_compound_log2fc_qvalue("dendrobatidis")
    assert result.empty


def test_load_compound_log2fc_qvalue_liq_vs_liq_stage_contrast_excluded(tmp_path, monkeypatch):
    fixture = tmp_path / "all_significant_features_summary.tsv"
    fixture.write_text(
        "row_id\tspecies\tcomparison\tlog2FC_a_over_b\tq_value\n"
        "12\tdendrobatidis\tdendrobatidis_liq_Zoospore_vs_liq_Developed\t8.0\t0.001\n"
    )
    monkeypatch.setattr(build_linkage_tables, "DIFFERENTIAL_FEATURES_PATH", fixture)
    result = build_linkage_tables.load_compound_log2fc_qvalue("dendrobatidis")
    assert result.empty


def test_strip_transcript_suffix():
    assert _strip_transcript_suffix("FCC698BD_000001-T1") == "FCC698BD_000001"
    assert _strip_transcript_suffix("F61BA062_000001-T2") == "F61BA062_000001"
    assert _strip_transcript_suffix("F61BA062_000001-mRNA-1") == "F61BA062_000001"
    assert _strip_transcript_suffix("FCC698BD_000001") == "FCC698BD_000001"


def test_compute_bgc_proteins_strips_transcript_suffix_before_coord_lookup():
    # pfam_calls protein_ids carry a "-T1" transcript suffix (real BFD PFAM
    # domtblout target-name format, e.g. "FCC698BD_000001-T1"); protein_coords
    # is keyed WITHOUT that suffix (real BFD .gbk CDS /locus_tag format, e.g.
    # "FCC698BD_000001"). Without suffix-stripping, every lookup misses and
    # has_bgc_context is False for every protein -- see
    # GENOME_BIOACTIVITY_LINKAGE.md C3.
    pfam_calls = pd.DataFrame(
        [{"protein_id": "FCC698BD_000001-T1", "pfam_id": "PF01397", "family": "terpene_synthase"}]
    )
    protein_coords = {"FCC698BD_000001": ("contig1", 1200, 2200)}
    regions = [{"record_id": "contig1", "start": 1000, "end": 5000}]

    bgc_proteins = compute_bgc_proteins(pfam_calls, protein_coords, regions)

    assert bgc_proteins == {"FCC698BD_000001-T1"}


def test_compute_bgc_proteins_raises_on_low_coord_resolution_rate():
    # If more than half of pfam_calls' unique protein_ids fail to resolve
    # against protein_coords even after suffix-stripping, this looks like a
    # namespace mismatch, not real biology -- fail loudly instead of
    # silently returning an empty/near-empty bgc_proteins set.
    pfam_calls = pd.DataFrame(
        [
            {"protein_id": "unknown1-T1", "pfam_id": "PF01397", "family": "terpene_synthase"},
            {"protein_id": "unknown2-T1", "pfam_id": "PF00109", "family": "pks"},
            {"protein_id": "unknown3-T1", "pfam_id": "PF00067", "family": "p450"},
        ]
    )
    protein_coords = {"unknown1": ("contig1", 1200, 2200)}
    regions = [{"record_id": "contig1", "start": 1000, "end": 5000}]

    try:
        compute_bgc_proteins(pfam_calls, protein_coords, regions)
        assert False, "expected ValueError"
    except ValueError as e:
        assert "namespace mismatch" in str(e) or "resolved" in str(e)
