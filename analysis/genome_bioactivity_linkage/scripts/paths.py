"""Shared path/locustag constants for the genome-bioactivity-linkage pipeline."""
from pathlib import Path

BFD_ROOT = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
REPO_ROOT = Path(__file__).resolve().parents[3]
GBL_ROOT = Path(__file__).resolve().parents[1]

SPECIES = {
    "dendrobatidis": {
        "out": "Batrachochytrium_dendrobatidis_JEL423",
        "locustag": "FCC698BD",
        "ncbi_assembly": "GCA_048537975.1_CMM_BatrDend_JEL423_V3",
    },
    "salamandrivorans": {
        "out": "Batrachochytrium_salamandrivorans_AMFP13",
        "locustag": "F61BA062",
        "ncbi_assembly": "GCA_002006685.2_Batr_sala_V2",
    },
}


def bfd_proteins_fasta(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    return BFD_ROOT / "genome_annotation" / out / "predict_results" / f"{out}.proteins.fa"


def bfd_gbk(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    return BFD_ROOT / "genome_annotation" / out / "predict_results" / f"{out}.gbk"


def bfd_antismash_json(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    base = BFD_ROOT / "genome_annotation" / out / "antismash_local"
    for name in (f"{out}.json", f"{out}.json.gz"):
        p = base / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"antiSMASH JSON not found for {species_key} under {base} "
        f"(expected {out}.json or {out}.json.gz) — BFD antiSMASH run may not be finished yet"
    )


def find_bfd_output(kind: str, species_key: str, suffix: str | None = None) -> Path:
    """Locate a locustag-bucketed BFD functional-annotation file, e.g. kind='pfam_hmmscan'.

    A locustag bucket for a given ``kind`` can legitimately contain more
    than one file for the same locustag -- e.g. every real ``pfam_hmmscan``
    bucket has BOTH ``<TAG>.domtblout.gz`` AND ``<TAG>.tblout.gz``. Without
    ``suffix``, ties are broken by ``sorted(matches)[0]`` (today this
    happens to pick the right file for pfam_hmmscan since 'd' < 't', but
    that's luck, not a guarantee, and picking the wrong file would silently
    produce a garbage/empty result with no error). Pass ``suffix`` (e.g.
    ``".domtblout.gz"``) to filter to only matches ending with it first.
    """
    locustag = SPECIES[species_key]["locustag"]
    search_root = BFD_ROOT / "results" / "function" / kind
    matches = sorted(search_root.glob(f"**/{locustag}*"))
    if suffix is not None:
        suffix_matches = [m for m in matches if str(m).endswith(suffix)]
        if suffix_matches:
            matches = suffix_matches
    if not matches:
        raise FileNotFoundError(
            f"No BFD '{kind}' output found for locustag {locustag} under {search_root} "
            f"— BFD functional-annotation run may not be finished yet"
        )
    return matches[0]
