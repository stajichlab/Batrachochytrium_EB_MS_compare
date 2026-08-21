"""Pfam-id -> biosynthetic domain family, and compound-class -> family mapping.

PF00494 (squalene/phytoene synthase) is deliberately its own family, not
folded into terpene_synthase -- see spec Stage 1.
"""

DOMAIN_FAMILIES: dict[str, set[str]] = {
    "terpene_synthase": {"PF01397", "PF03936"},
    "squalene_phytoene_synthase": {"PF00494"},
    "pks": {"PF00109", "PF08659", "PF00698"},
    "nrps": {"PF00668", "PF00501", "PF00550"},
    "dmats_prenyltransferase": {"PF11991"},
    "p450": {"PF00067"},
}

# Only compound classes with a defined candidate domain family are mapped;
# unmapped classes (e.g. "Fatty acyls") are deliberately absent so Stage 4
# leaves them out of the candidate table rather than force-fitting them.
COMPOUND_CLASS_TO_FAMILY: dict[str, str] = {
    "Terpenoids": "terpene_synthase",
    "Polyketides": "pks",
    "Alkaloids (linear polyketides)": "pks",
    "Amino acids and Peptides": "nrps",
}

_PFAM_TO_FAMILY = {pfam: family for family, pfams in DOMAIN_FAMILIES.items() for pfam in pfams}


def classify_pfam(pfam_id: str) -> str | None:
    """Map a bare Pfam accession (e.g. 'PF01397', no version suffix) to its family, or None."""
    return _PFAM_TO_FAMILY.get(pfam_id)
