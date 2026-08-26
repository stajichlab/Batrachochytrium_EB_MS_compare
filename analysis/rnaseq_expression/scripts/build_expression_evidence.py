#!/usr/bin/env python3
"""Fold RNA-seq expression evidence into the genome-bioactivity candidate tables.

The STAR --quantMode GeneCounts step produced no usable counts (the raw NCBI
GFF3 exon lines carry no gene_id attribute), so reads-per-gene are recomputed
with featureCounts from a gffread-converted GTF (see run_featurecounts.sh,
RNASEQ_EXPRESSION.md). This script:

1. reads per-sample fragment counts (counts_s0.txt, the best-assigning
   strandedness mode for these public runs),
2. maps gene_id -> NCBI locus_tag (strip the "gene-" prefix),
3. maps each candidate BFD protein -> its reciprocal-best-hit NCBI protein
   -> that protein's locus_tag,
4. joins the candidate table to per-gene expression, adding
   rna_total_raw (gene-level raw counts summed across that species' runs),
   rna_n_rep_cov (how many of the replicates have >=min_counts reads),
   rna_is_expressed (rna_n_rep_cov >= min_reps).

Caveat carried from the RNA-seq pipeline design: these public RNA-seq runs are
NOT from the liquid-culture growth condition sampled by the mass-spec data, so
this is presence/absence expression evidence, not condition-matched
co-expression. Absence is weaker evidence than presence because of that.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
GBL = REPO / "analysis" / "genome_bioactivity_linkage"
RSEQ = REPO / "analysis" / "rnaseq_expression"

SPECIES_SAMPLES = {
    "dendrobatidis": ["SRR27683881", "SRR27683880", "SRR27683879"],
    "salamandrivorans": ["SRR13012113", "SRR13012117", "SRR13012121", "SRR13012125", "SRR13012129"],
}

MIN_COUNTS = 10  # raw fragment-count floor for "detected" in a sample
MIN_REPS = 1  # how many replicates must meet MIN_COUNTS to call "expressed"


def load_featurecounts(species: str) -> pd.DataFrame:
    f = RSEQ / "results" / "gene_counts" / species / "counts_s0.txt"
    df = pd.read_csv(f, sep="\t", comment="#", skiprows=1)
    df = df.rename(columns={"Geneid": "gene_id"})
    df["locus_tag"] = df["gene_id"].str.replace(r"^gene-", "", regex=True)
    cols = ["gene_id", "locus_tag"] + [c for c in df.columns if c.endswith(".bam")]
    return df[cols]


def load_rbh(species: str) -> pd.DataFrame:
    rbh_dir = GBL / "results" / "rbh" / species
    fwd = pd.read_csv(rbh_dir / "fwd.tsv", sep="\t", names=["bfd", "ref", "pi", "len", "mm", "go", "qs", "qe", "ss", "se", "ev", "bs"])
    rev = pd.read_csv(rbh_dir / "rev.tsv", sep="\t", names=["ref", "bfd", "pi", "len", "mm", "go", "qs", "qe", "ss", "se", "ev", "bs"])
    fwd_f = fwd.sort_values("bs", ascending=False).drop_duplicates("bfd")
    rev_f = rev.sort_values("bs", ascending=False).drop_duplicates("ref")
    seen = set()
    rows = []
    for _, r in fwd_f.iterrows():
        if r.ref in rev_f.index and rev_f.loc[r.ref, "bfd"] == r.bfd and r.bfd not in seen:
            seen.add(r.bfd)
            rows.append({"bfd_protein_id": r.bfd, "reference_protein_id": r.ref})
    return pd.DataFrame(rows, columns=["bfd_protein_id", "reference_protein_id"])


def main() -> None:
    for species, samples in SPECIES_SAMPLES.items():
        counts = load_featurecounts(species)
        rbh = load_rbh(species)

        prot2loc = {}
        with open(REPO / "analysis" / "genome_bioactivity_linkage" / "results" / "reference_annotation" / species / "genomic.gff") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.rstrip().split("\t")
                if len(f) < 9 or f[2] != "CDS":
                    continue
                a = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
                pid, loc = a.get("protein_id"), a.get("locus_tag")
                if pid and loc:
                    prot2loc[pid] = loc

        counts["gene_total_raw"] = counts[[c for c in counts.columns if c.endswith(".bam")]].sum(axis=1)
        counts["n_rep_ge_min"] = (counts[[c for c in counts.columns if c.endswith(".bam")]] >= MIN_COUNTS).sum(axis=1)

        cand = pd.read_csv(GBL / "results" / f"{species}_candidate_table.tsv", sep="\t")
        # drop columns this script itself wrote on a previous pass so the
        # merge below doesn't suffix them into *_x/_y duplicates
        for drop in ["reference_protein_id", "ref_locus", "gene_total_raw",
                     "n_rep_ge_min", "rna_is_expressed", "rna_no_evidence"]:
            cand = cand.drop(columns=[drop], errors="ignore")
        # map candidate BFD protein -> ref protein (top-hit or RBH; prefer RBH)
        ref_map = dict(zip(rbh.bfd_protein_id, rbh.reference_protein_id))
        fwd_top = pd.read_csv(GBL / "results" / "rbh" / species / "fwd.tsv", sep="\t", names=["bfd", "ref", "pi", "len", "mm", "go", "qs", "qe", "ss", "se", "ev", "bs"]).sort_values("bs", ascending=False).drop_duplicates("bfd")
        top_map = dict(zip(fwd_top.bfd, fwd_top.ref))

        cand["reference_protein_id"] = cand["candidate_protein_id"].map(lambda p: ref_map.get(p) or top_map.get(p))
        cand["ref_locus"] = cand["reference_protein_id"].map(prot2loc)
        expr = counts.set_index("locus_tag")[["gene_total_raw", "n_rep_ge_min"]]
        cand = cand.merge(expr, left_on="ref_locus", right_index=True, how="left")
        cand["rna_is_expressed"] = cand["n_rep_ge_min"].fillna(0).astype(int) >= MIN_REPS
        cand["rna_is_expressed"] = cand["rna_is_expressed"].astype(bool)
        # expression flag missing entirely when no RBH/top hit maps to a locus
        cand["rna_no_evidence"] = cand["ref_locus"].isna()

        out_cand = GBL / "results" / f"{species}_candidate_table.tsv"
        cand.to_csv(out_cand, sep="\t", index=False)

        # per-gene expression summary for downstream use
        counts = counts.merge(rbh, left_on="locus_tag", right_on="reference_protein_id", how="left")
        rbh_expr = counts.groupby("locus_tag", as_index=False)[["gene_total_raw", "n_rep_ge_min"]].max()
        expr_out = RSEQ / "results" / "gene_counts" / species / "gene_expression.tsv"
        counts.to_csv(expr_out, sep="\t", index=False)

        total = len(cand)
        n_expr = int(cand["rna_is_expressed"].sum())
        n_expr_true = int((cand["rna_is_expressed"] & (~cand["is_extracellular"])).sum())
        print(
            f"{species}: candidates={total} rna_expressed={n_expr} "
            f"(extracellular={n_expr_true}, no_locus_mapping={int(cand['rna_no_evidence'].sum())}) "
            f"-> {out_cand}"
        )


if __name__ == "__main__":
    main()
