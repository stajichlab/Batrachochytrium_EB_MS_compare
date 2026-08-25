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
    """Locate the antiSMASH JSON for a species.

    Prefers BFD's own antiSMASH sub-run (never materialized as of
    2026-08-25 for either Batrachochytrium genome). Falls back to this
    project's own antiSMASH run on the NCBI reference GBFF
    (run_antismash_reference.sh) — confirmed safe to use for
    `has_bgc_context` because the NCBI assembly and BFD's own gene-model
    assembly share identical contig ids/lengths (e.g. Bd `CP161923.1` =
    4,539,083 bp, Bsal `LYON02000001.1` = 5,593,528 bp in both; 20/165
    total records match exactly), so BGC-region coordinates from this run
    are directly comparable to BFD's own gene coordinates.
    """
    out = SPECIES[species_key]["out"]
    base = BFD_ROOT / "genome_annotation" / out / "antismash_local"
    for name in (f"{out}.json", f"{out}.json.gz"):
        p = base / name
        if p.exists():
            return p

    fallback = GBL_ROOT / "results" / "antismash_ncbi" / species_key / f"{out}.json"
    if fallback.exists():
        return fallback

    raise FileNotFoundError(
        f"antiSMASH JSON not found for {species_key} under {base} "
        f"(expected {out}.json or {out}.json.gz) — BFD antiSMASH run may not be finished yet — "
        f"nor under the fallback {fallback} (run run_antismash_reference.sh to produce it)"
    )


def find_bfd_output(kind: str, species_key: str, suffix: str | None = None) -> Path:
    """Locate a locustag-bucketed functional-annotation file, e.g. kind='pfam_hmmscan'.

    Prefers BFD's own shared functional-annotation run
    (``BFD_ROOT/results/function/<kind>/``). Falls back to this project's
    own local run (``GBL_ROOT/results/<kind>/<species_key>/``, produced by
    ``run_pfam.sh`` / ``run_signalp.sh`` / ``run_predgpi.sh``) when BFD's
    shared run has not reached this locustag yet — as of 2026-08-25, BFD's
    ``pfam_hmmscan``/``signalp``/``predgpi`` have no output at all for
    either Batrachochytrium genome (confirmed: zero files matching either
    locustag under BFD_ROOT). The local run reuses BFD's own tool
    invocations (same Pfam-A.hmm database, same ``hmmsearch --cut_ga``,
    same SignalP6/PredGPI CLI flags) and writes the identical filename/
    format conventions, so the parsers in ``parse_pfam_domains.py`` /
    ``merge_secretion.py`` need no changes regardless of which source
    supplied the file.

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
    search_roots = [
        BFD_ROOT / "results" / "function" / kind,
        GBL_ROOT / "results" / kind / species_key,
    ]
    matches: list[Path] = []
    for search_root in search_roots:
        matches = sorted(search_root.glob(f"**/{locustag}*"))
        if matches:
            break
    if not matches:
        raise FileNotFoundError(
            f"No '{kind}' output found for locustag {locustag} under any of "
            f"{[str(r) for r in search_roots]} — neither BFD's shared run nor "
            f"this project's local run (run_pfam.sh/run_signalp.sh/run_predgpi.sh) "
            f"has produced it yet"
        )
    if suffix is not None:
        suffix_matches = [m for m in matches if str(m).endswith(suffix)]
        if not suffix_matches:
            raise FileNotFoundError(
                f"No '{kind}' output for locustag {locustag} under {search_root} "
                f"matches suffix {suffix!r} (found: {[m.name for m in matches]})"
            )
        matches = suffix_matches
    return matches[0]
