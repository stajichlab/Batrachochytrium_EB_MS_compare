#!/usr/bin/env python3
"""Secreted-protease candidate table (companion to build_linkage_tables.py's
compound-linked biosynthetic-gene candidates).

Proteases don't biosynthesize a small-molecule compound the way
NRPS/PKS/terpene-synthase genes do, so this is NOT joined to SIRIUS
compound classes the way the main candidate table is -- there is no
"protease produces this metabolite" relationship to encode. Instead this
ranks BFD-predicted proteins with a MEROPS peptidase-family hit by the
same secretion/orthology/expression evidence used elsewhere in this
pipeline, on the premise that a SECRETED, expressed, orthology-confirmed
protease is a plausible extracellular virulence factor (fungal pathogens,
including chytrids, commonly use secreted proteases to degrade host
proteins) -- and a candidate explanation for some of the "Amino acids and
Peptides"-class liquid-enriched metabolomics signal (proteolytic
fragments of host/media proteins), independent of the NRPS hypothesis in
TIER1_NRPS_CHARACTERIZATION.md.

Reuses this pipeline's existing secretion (SignalP+DeepTMHMM+PredGPI),
RBH ortholog, and RNA-seq expression machinery verbatim -- see
build_linkage_tables.run_for_species / build_expression_evidence.py for
the same joins on the compound-linkage side.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_secretion import load_predgpi_gff3, load_signalp_gff3, predicted_extracellular
from parse_deeptmhmm import parse_tmrs_gff3
from parse_merops import best_merops_hit, load_merops_blasttab, load_merops_families
from parse_rbh import reciprocal_best_hits
from paths import GBL_ROOT, SPECIES, find_bfd_output

MIN_COUNTS = 10  # matches build_expression_evidence.py's raw-count floor
MIN_REPS = 1

SPECIES_SAMPLES = {
    "dendrobatidis": ["SRR27683881", "SRR27683880", "SRR27683879"],
    "salamandrivorans": ["SRR13012113", "SRR13012117", "SRR13012121", "SRR13012125", "SRR13012129"],
}


def load_expression(species_key: str, rbh: pd.DataFrame) -> pd.DataFrame:
    """Per-BFD-protein expression evidence via the RBH -> NCBI locus ->
    featureCounts join (same approach as build_expression_evidence.py)."""
    repo = GBL_ROOT.parents[1]
    counts_path = repo / "analysis" / "rnaseq_expression" / "results" / "gene_counts" / species_key / "counts_s0.txt"
    gff_path = GBL_ROOT / "results" / "reference_annotation" / species_key / "genomic.gff"
    if not counts_path.exists() or not gff_path.exists():
        return pd.DataFrame(columns=["protein_id", "reference_protein_id", "ref_locus", "gene_total_raw", "n_rep_ge_min", "rna_is_expressed"])

    counts = pd.read_csv(counts_path, sep="\t", comment="#", skiprows=1).rename(columns={"Geneid": "gene_id"})
    counts["locus_tag"] = counts["gene_id"].str.replace(r"^gene-", "", regex=True)
    bam_cols = [c for c in counts.columns if c.endswith(".bam")]
    counts["gene_total_raw"] = counts[bam_cols].sum(axis=1)
    counts["n_rep_ge_min"] = (counts[bam_cols] >= MIN_COUNTS).sum(axis=1)

    prot2loc: dict[str, str] = {}
    with open(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip().split("\t")
            if len(f) < 9 or f[2] != "CDS":
                continue
            attrs = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
            pid, loc = attrs.get("protein_id"), attrs.get("locus_tag")
            if pid and loc:
                prot2loc[pid] = loc

    expr = pd.DataFrame({"reference_protein_id": rbh["reference_protein_id"], "protein_id": rbh["bfd_protein_id"]})
    expr["ref_locus"] = expr["reference_protein_id"].map(prot2loc)
    counts_by_locus = counts.set_index("locus_tag")[["gene_total_raw", "n_rep_ge_min"]]
    expr = expr.merge(counts_by_locus, left_on="ref_locus", right_index=True, how="left")
    expr["rna_is_expressed"] = (expr["n_rep_ge_min"].fillna(0).astype(int) >= MIN_REPS)
    return expr


def run_for_species(species_key: str) -> pd.DataFrame:
    merops_path = find_bfd_output("merops", species_key)
    families = load_merops_families()
    merops = best_merops_hit(load_merops_blasttab(merops_path), families)
    merops = merops.rename(columns={"family": "merops_family", "clan": "merops_clan"})
    # Peptidase-INHIBITOR family hits ("I" catalytic type) are not
    # peptidases themselves -- flag rather than silently drop, so the
    # table stays auditable (same "informational column, not a hard
    # filter" precedent as require_extracellular=False elsewhere in this
    # pipeline).
    merops["is_inhibitor_family"] = merops["catalytic_type"] == "inhibitor"

    deeptmhmm_dir = GBL_ROOT / "results" / "deeptmhmm" / species_key
    deeptmhmm_gff3 = parse_tmrs_gff3(deeptmhmm_dir / "TMRs.gff3")
    signalp_path = find_bfd_output("signalp", species_key, suffix=".signalp.gff3.gz")
    predgpi_path = find_bfd_output("predgpi", species_key)
    signalp = load_signalp_gff3(signalp_path)
    predgpi = load_predgpi_gff3(predgpi_path)
    secretion = predicted_extracellular(signalp, deeptmhmm_gff3, predgpi)

    rbh_dir = GBL_ROOT / "results" / "rbh" / species_key
    rbh = reciprocal_best_hits(rbh_dir / "fwd.tsv", rbh_dir / "rev.tsv")
    rbh_confirmed = set(rbh["bfd_protein_id"])

    table = merops.merge(secretion, on="protein_id", how="left")
    table["is_extracellular"] = table["is_extracellular"].fillna(False).astype(bool)
    table["is_cross_ref_confirmed"] = table["protein_id"].isin(rbh_confirmed)

    expr = load_expression(species_key, rbh)
    table = table.merge(expr.drop(columns=["reference_protein_id"], errors="ignore"), on="protein_id", how="left")
    table["rna_is_expressed"] = table["rna_is_expressed"].fillna(False).astype(bool)

    table = table.sort_values(
        ["is_extracellular", "rna_is_expressed", "evalue"], ascending=[False, False, True]
    ).reset_index(drop=True)
    return table


if __name__ == "__main__":
    out_dir = GBL_ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    for species_key in SPECIES:
        table = run_for_species(species_key)
        out_path = out_dir / f"{species_key}_protease_candidates.tsv"
        table.to_csv(out_path, sep="\t", index=False)
        n_secreted = int(table["is_extracellular"].sum())
        n_secreted_expr = int((table["is_extracellular"] & table["rna_is_expressed"]).sum())
        print(
            f"{species_key}: {len(table)} MEROPS-hit proteins, {n_secreted} secreted-candidate "
            f"({n_secreted_expr} also expressed) -> {out_path}"
        )
