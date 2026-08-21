import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import domain_families  # noqa: E402


def test_terpene_synthase_excludes_squalene_phytoene_synthase():
    assert domain_families.classify_pfam("PF01397") == "terpene_synthase"
    assert domain_families.classify_pfam("PF03936") == "terpene_synthase"
    assert domain_families.classify_pfam("PF00494") == "squalene_phytoene_synthase"
    assert domain_families.classify_pfam("PF00494") != "terpene_synthase"


def test_pks_and_nrps_and_unknown():
    assert domain_families.classify_pfam("PF00109") == "pks"
    assert domain_families.classify_pfam("PF00668") == "nrps"
    assert domain_families.classify_pfam("PF99999") is None


def test_compound_class_to_family_map_omits_unmapped_classes():
    assert domain_families.COMPOUND_CLASS_TO_FAMILY["Terpenoids"] == "terpene_synthase"
    assert domain_families.COMPOUND_CLASS_TO_FAMILY["Polyketides"] == "pks"
    assert "Fatty acyls" not in domain_families.COMPOUND_CLASS_TO_FAMILY
