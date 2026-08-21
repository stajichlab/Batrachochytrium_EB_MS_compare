# analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py
"""End-to-end driver: BFD Stage-1/2 outputs + this repo's Stage-3/4 modules
-> per-species tiered candidate tables. See GENOME_BIOACTIVITY_LINKAGE.md.

Two deliberate deviations from a naive first-pass implementation (see
GENOME_BIOACTIVITY_LINKAGE.md "Known caveats" / task-10 review notes):

1. ``has_bgc_context`` is NOT computed via fullhmmer-hit membership --
   ``full_hmmer`` is antiSMASH's genome-wide internal PFAM scan, not a
   BGC-region membership signal (most antiSMASH JSON records have an empty
   ``areas`` list yet still contribute full_hmmer hits). Instead each PFAM
   call's protein is looked up by its REAL genomic coordinates (parsed once
   from the BFD ``.gbk`` file) and tested against antiSMASH's BGC regions via
   ``parse_antismash_clusters.protein_in_bgc``.
2. ``compounds["log2fc"]``/``["q_value"]`` are NOT hardcoded placeholders --
   they are joined from ``analysis/differential_features_primary/
   all_significant_features_summary.tsv`` (liq-fraction contrasts only,
   lowest-q-value row kept per compound). Compounds with no matching row in
   that file are excluded from the candidate table entirely.
"""
import sys
from pathlib import Path

import pandas as pd
from Bio import SeqIO

sys.path.insert(0, str(Path(__file__).resolve().parent))
from background_subtraction import fungal_over_blank_ratio, load_feature_intensities, load_metadata
from domain_families import COMPOUND_CLASS_TO_FAMILY
from link_compounds_to_genes import build_candidate_table
from merge_secretion import predicted_extracellular
from parse_antismash_clusters import load_fullhmmer_hits, load_regions, protein_in_bgc
from parse_deeptmhmm import parse_tmrs_gff3
from parse_pfam_domains import classify_domains, parse_domtblout
from parse_rbh import reciprocal_best_hits
from paths import GBL_ROOT, SPECIES, bfd_antismash_json, bfd_gbk, find_bfd_output

DIFFERENTIAL_FEATURES_PATH = (
    GBL_ROOT.parents[0] / "differential_features_primary" / "all_significant_features_summary.tsv"
)


def build_gene_domain_table(
    pfam_calls: pd.DataFrame,
    bgc_proteins: set,
    secretion: pd.DataFrame,
    rbh_confirmed_proteins: set,
) -> pd.DataFrame:
    merged = pfam_calls.merge(secretion, on="protein_id", how="left")
    merged["has_bgc_context"] = merged["protein_id"].isin(bgc_proteins)
    merged["is_cross_ref_confirmed"] = merged["protein_id"].isin(rbh_confirmed_proteins)
    merged["is_extracellular"] = merged["is_extracellular"].fillna(False)
    return merged[
        [
            "protein_id", "family", "has_bgc_context", "is_cross_ref_confirmed", "is_extracellular",
        ]
    ]


def load_protein_coords(gbk_path: Path) -> dict:
    """Parse a BFD .gbk file once into {locus_tag: (record_id, start, end)}.

    ``feature.location.start``/``.end`` are 0-based genomic coordinates on the
    parent contig record -- the same convention antiSMASH's own region
    coordinates use (verified directly against the antiSMASH JSON: record
    ``.id`` values match GenBank record ids, e.g. ``CP161923.1``). These are
    NOT the same as full_hmmer's ``protein_start``/``protein_end``, which are
    amino-acid offsets within the protein sequence.
    """
    coords = {}
    for record in SeqIO.parse(str(gbk_path), "genbank"):
        for feature in record.features:
            if feature.type != "CDS":
                continue
            locus_tags = feature.qualifiers.get("locus_tag")
            if not locus_tags:
                continue
            locus_tag = locus_tags[0]
            coords[locus_tag] = (record.id, int(feature.location.start), int(feature.location.end))
    return coords


def compute_bgc_proteins(pfam_calls: pd.DataFrame, protein_coords: dict, regions: list) -> set:
    """Per-protein BGC-region membership via real genomic coordinates.

    Deliberately NOT ``{p for p in fullhmmer["protein_id"] if p in
    set(pfam_calls["protein_id"])}`` -- see module docstring deviation 1.
    """
    bgc_proteins = set()
    for protein_id in pfam_calls["protein_id"].unique():
        coord = protein_coords.get(protein_id)
        if coord is None:
            continue
        record_id, start, end = coord
        if protein_in_bgc(protein_id, (start, end), regions, record_id):
            bgc_proteins.add(protein_id)
    return bgc_proteins


def load_compound_log2fc_qvalue(species_key: str) -> pd.DataFrame:
    """Join compound row_ids to their most-significant liq-fraction contrast
    result (lowest q_value) from the primary differential-features table.

    Deliberately NOT a hardcoded ``log2fc=0.0``/``q_value=1.0`` placeholder --
    see module docstring deviation 2. A compound with no matching liq-fraction
    contrast row is excluded upstream (never significant in the liquid
    fraction, so it should not receive a fabricated tie-breaker value).
    """
    diff = pd.read_csv(
        DIFFERENTIAL_FEATURES_PATH,
        sep="\t",
        usecols=["row_id", "species", "comparison", "log2FC_a_over_b", "q_value"],
    )
    diff = diff[(diff["species"] == species_key) & (diff["comparison"].str.contains("liq"))]
    diff = diff.sort_values("q_value").drop_duplicates("row_id", keep="first")
    return diff.rename(columns={"log2FC_a_over_b": "log2fc"})[["row_id", "log2fc", "q_value"]]


def run_for_species(species_key: str) -> pd.DataFrame:
    domtbl_path = find_bfd_output("pfam_hmmscan", species_key)
    pfam_calls = classify_domains(parse_domtblout(domtbl_path))

    antismash_json = bfd_antismash_json(species_key)
    regions = load_regions(antismash_json)
    load_fullhmmer_hits(antismash_json)  # not used for has_bgc_context -- see deviation 1
    protein_coords = load_protein_coords(bfd_gbk(species_key))
    bgc_proteins = compute_bgc_proteins(pfam_calls, protein_coords, regions)

    rbh_dir = GBL_ROOT / "results" / "rbh" / species_key
    rbh = reciprocal_best_hits(rbh_dir / "fwd.tsv", rbh_dir / "rev.tsv")
    rbh_confirmed_proteins = set(rbh["bfd_protein_id"])

    deeptmhmm_dir = GBL_ROOT / "results" / "deeptmhmm" / species_key
    deeptmhmm_gff3 = parse_tmrs_gff3(deeptmhmm_dir / "TMRs.gff3")
    signalp_path = find_bfd_output("signalp", species_key)
    predgpi_path = find_bfd_output("predgpi", species_key)
    signalp = pd.read_csv(signalp_path, sep="\t")
    predgpi = pd.read_csv(predgpi_path, sep="\t")
    secretion = predicted_extracellular(signalp, deeptmhmm_gff3, predgpi)

    gene_domains = build_gene_domain_table(pfam_calls, bgc_proteins, secretion, rbh_confirmed_proteins)

    meta = load_metadata()
    features = load_feature_intensities()
    sirius = pd.read_csv(GBL_ROOT.parents[0] / "sirius_annotation" / "sirius_annotations.tsv", sep="\t")

    species_full_name = f"Batrachochytrium {species_key}"
    life_stages = ["Zoospore", "Sporangium", "Mature"]
    passing_ids = set()
    for stage in life_stages:
        ratio = fungal_over_blank_ratio(features, meta, species=species_full_name, life_stage=stage, min_fc=2.0)
        passing_ids |= set(ratio.loc[ratio["passes_background_filter"], "row_id"])

    sirius = sirius[sirius["row ID"].isin(passing_ids)].copy()
    sirius["compound_class"] = sirius["sirius_npc_pathway"].where(
        sirius["sirius_npc_pathway"].isin(COMPOUND_CLASS_TO_FAMILY), sirius["sirius_npc_class"]
    )
    compounds = sirius.rename(columns={"row ID": "row_id"})[["row_id", "compound_class"]].copy()

    log2fc_qvalue = load_compound_log2fc_qvalue(species_key)
    compounds = compounds.merge(log2fc_qvalue, on="row_id", how="inner")

    return build_candidate_table(compounds, gene_domains)


if __name__ == "__main__":
    out_dir = GBL_ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    for species_key in SPECIES:
        table = run_for_species(species_key)
        out_path = out_dir / f"{species_key}_candidate_table.tsv"
        table.to_csv(out_path, sep="\t", index=False)
        print(f"{species_key}: {len(table)} candidate rows -> {out_path}")
