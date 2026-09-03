#!/usr/bin/env python3
"""Molecular-family (component) level analysis of the FBMN network.

Why this tier exists
--------------------
GOALS.md goals 1 and 4 (characterize molecular families; probe DeltaMZ
ladders) were never attempted. They matter more after the 2026-09-02
corrections than before, because per-feature evidence in this dataset is
weak by construction: n=5 per group, a BH threshold that only fires when
~16% of features separate at once, and a medium background that mimics the
biological hypothesis.

A molecular family is a different and stronger unit of evidence. If a
homologous series -- features linked by MS2 cosine, differing by a
chemically interpretable DeltaMZ -- is *entirely* blank-clearing, that is
much harder to explain as medium background than any single feature is,
because the medium would have to contribute the whole series.

What it does
------------
1. Builds connected components from `filtered_pairs.tsv` (MS2 cosine edges).
2. Restricts to features that survived the artifact filter AND are present
   in the analysis matrix, so components are described in the same feature
   universe every other tier uses. Component sizes are reported both before
   and after this restriction, since dropping isotope/adduct rows is
   expected to shrink families substantially -- those rows are exactly what
   inflates a "molecular family" with non-independent copies of one molecule.
3. Annotates each component with: size, per-species blank-clearing member
   counts, SIRIUS class composition, and whether it is fully blank-clearing.
4. Tests, per component, whether its blank-clearing fraction exceeds the
   genome-wide background rate (Fisher exact, BH-corrected). A component
   that is significantly blank-enriched is a candidate fungal-derived family.
5. Classifies edge DeltaMZ against a curated table of common homologous
   steps (CH2, H2, H2O, C2H4, O, NH, CH3, CO, C2H2) within tolerance, so
   ladders can be read off directly.

Interpretation limits
---------------------
- Components come from the GNPS2 run's own cosine threshold; this script does
  not re-network anything.
- Edges connect features, not necessarily true analogs: adduct/in-source
  pairs also cosine-match. The artifact filter removes the systematic classes
  (isotopes, non-default adducts, ISF) but co-isolation chimeras remain.
- `all_blank_clearing` is evidence about the FAMILY, not proof any member is
  a secreted secondary metabolite. The peptide-origin caveat still applies:
  a family of medium-protein digest peptides can be entirely blank-clearing
  if the fungus generates all of them.

Usage: python3 scripts/component_analysis.py
Outputs: analysis/molecular_network/{components.tsv, delta_mz_ladders.tsv,
         component_summary.md, component_sizes.png}
"""
from __future__ import annotations

import sys
from collections import defaultdict, Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO = Path(__file__).resolve().parents[3]
PAIRS = REPO / "data" / "raw" / "gnps2_e9838293_bagel" / "nf_output" / "networking" / "filtered_pairs.tsv"
LINKED = REPO / "analysis" / "ordination" / "linked_data"
SIRIUS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT = REPO / "analysis" / "molecular_network"

sys.path.insert(0, str(REPO / "analysis" / "genome_bioactivity_linkage" / "scripts"))
from background_subtraction import (  # noqa: E402
    fungal_over_blank_ratio,
    load_feature_intensities,
    load_metadata,
)

SPECIES = ["Batrachochytrium dendrobatidis", "Batrachochytrium salamandrivorans"]
STAGES = ["Zoospore", "Sporangium", "Mature"]

# Common homologous / modification steps. Tolerance is the workflow's own
# fragment tolerance (0.05 Da, see GOALS.md workflow params).
DELTA_STEPS = {
    "CH2 (homolog)": 14.0157,
    "2xCH2": 28.0313,
    "H2 (saturation)": 2.0157,
    "H2O": 18.0106,
    "C2H4": 28.0313,
    "O (oxidation)": 15.9949,
    "NH": 15.0109,
    "CH3": 15.0235,
    "CO": 27.9949,
    "C2H2": 26.0157,
}
DELTA_TOL = 0.05


def components_from_edges(edges: pd.DataFrame) -> dict[int, int]:
    """feature_id -> component_id, via union over cosine edges."""
    adj = defaultdict(set)
    for a, b in zip(edges["CLUSTERID1"], edges["CLUSTERID2"]):
        adj[a].add(b)
        adj[b].add(a)
    comp: dict[int, int] = {}
    cid = 0
    for node in adj:
        if node in comp:
            continue
        stack, members = [node], []
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp[x] = cid
            members.append(x)
            stack.extend(y for y in adj[x] if y not in comp)
        cid += 1
    return comp


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if not PAIRS.exists():
        sys.exit(f"missing {PAIRS}")
    edges = pd.read_csv(PAIRS, sep="\t")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz", usecols=["row_id", "mz", "rt"])
    kept = set(feat["row_id"].astype(int))
    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    if meta["is_C_companion"].any():
        sys.exit("media blanks present in sample_metadata.csv -- rebuild linked_data")

    comp = components_from_edges(edges)
    n_nodes_all = len(comp)
    n_comp_all = len(set(comp.values()))
    print(f"network: {len(edges)} edges, {n_nodes_all} nodes, {n_comp_all} components",
          file=sys.stderr)
    # Cross-check our union-find against GNPS2's own ComponentIndex: every edge
    # must join two nodes we placed in the SAME component, and the partition
    # must agree 1:1. A mismatch means our traversal is wrong, not GNPS2's.
    if "ComponentIndex" in edges.columns:
        bad = sum(1 for a, b in zip(edges["CLUSTERID1"], edges["CLUSTERID2"])
                  if comp[a] != comp[b])
        pairing = edges.assign(ours=[comp[a] for a in edges["CLUSTERID1"]])
        n_ours_per_gnps = pairing.groupby("ComponentIndex")["ours"].nunique()
        n_gnps_per_ours = pairing.groupby("ours")["ComponentIndex"].nunique()
        print(f"  cross-check vs GNPS2 ComponentIndex: {bad} edges spanning our "
              f"components (must be 0); max our-components per GNPS2 component "
              f"{n_ours_per_gnps.max()}, max GNPS2 per ours {n_gnps_per_ours.max()} "
              f"(both must be 1)", file=sys.stderr)
        if bad or n_ours_per_gnps.max() != 1 or n_gnps_per_ours.max() != 1:
            sys.exit("component traversal disagrees with GNPS2 ComponentIndex")

    # Blank-clearing sets, per species (paired-by-plate rule, union over stages).
    features_raw = load_feature_intensities()
    meta_raw = load_metadata()
    blank_ok = {}
    for sp in SPECIES:
        ids: set[int] = set()
        for st in STAGES:
            r = fungal_over_blank_ratio(features_raw, meta_raw, species=sp, life_stage=st)
            ids |= {int(x) for x in r.loc[r["passes_background_filter"], "row_id"]}
        blank_ok[sp.split()[-1]] = ids
        print(f"  blank-clearing {sp.split()[-1]}: {len(ids)}", file=sys.stderr)

    sirius = pd.read_csv(SIRIUS, sep="\t")
    sirius["row ID"] = sirius["row ID"].astype(int)
    npc = dict(zip(sirius["row ID"], sirius["sirius_npc_pathway"]))
    struct = dict(zip(sirius["row ID"], sirius["sirius_structure_name"]))

    # Background rate = fraction of ANALYSIS-MATRIX features that clear the blank.
    bg = {s: len(v & kept) / max(len(kept), 1) for s, v in blank_ok.items()}
    print(f"  background blank-clearing rate: "
          f"{ {k: round(v,4) for k,v in bg.items()} }", file=sys.stderr)

    members = defaultdict(list)
    for node, cid in comp.items():
        members[cid].append(int(node))

    rows = []
    for cid, mem in members.items():
        in_matrix = [m for m in mem if m in kept]
        if not in_matrix:
            continue
        rec = {
            "component": cid,
            "size_network": len(mem),
            "size_in_analysis_matrix": len(in_matrix),
            "n_sirius_annotated": sum(1 for m in in_matrix if m in npc and pd.notna(npc[m])),
        }
        classes = Counter(npc[m] for m in in_matrix if m in npc and pd.notna(npc[m]))
        rec["dominant_npc_pathway"] = classes.most_common(1)[0][0] if classes else ""
        names = [str(struct[m]) for m in in_matrix if m in struct and pd.notna(struct[m])]
        rec["example_structures"] = "; ".join(names[:3])
        for sp_short, ids in blank_ok.items():
            k = sum(1 for m in in_matrix if m in ids)
            n = len(in_matrix)
            rec[f"n_blank_clearing_{sp_short}"] = k
            rec[f"frac_blank_clearing_{sp_short}"] = round(k / n, 4)
            rec[f"all_blank_clearing_{sp_short}"] = (k == n and n >= 3)
            # Fisher: this component's blank-clearing members vs the background rate.
            bg_k = len(ids & kept)
            bg_n = len(kept)
            if n >= 3:
                _, p = fisher_exact([[k, n - k], [bg_k - k, (bg_n - n) - (bg_k - k)]],
                                    alternative="greater")
            else:
                p = 1.0
            rec[f"fisher_p_{sp_short}"] = p
        rows.append(rec)

    df = pd.DataFrame(rows).sort_values("size_in_analysis_matrix", ascending=False)
    for sp_short in blank_ok:
        p = df[f"fisher_p_{sp_short}"].to_numpy()
        order = np.argsort(p)
        q = np.empty_like(p)
        ranked = p[order] * len(p) / (np.arange(len(p)) + 1)
        ranked = np.minimum.accumulate(ranked[::-1])[::-1]
        q[order] = np.clip(ranked, 0, 1)
        df[f"q_blank_enriched_{sp_short}"] = q
    df.to_csv(OUT / "components.tsv", sep="\t", index=False)

    # DeltaMZ ladder classification on edges whose BOTH ends survive filtering.
    e = edges[edges["CLUSTERID1"].isin(kept) & edges["CLUSTERID2"].isin(kept)].copy()
    e["abs_delta"] = e["DeltaMZ"].abs()
    def classify(d):
        for label, target in DELTA_STEPS.items():
            if abs(d - target) <= DELTA_TOL:
                return label
        return ""
    e["delta_step"] = e["abs_delta"].map(classify)
    lad = (e[e["delta_step"] != ""]
           .groupby("delta_step")
           .agg(n_edges=("delta_step", "size"),
                median_cosine=("Cosine", "median"))
           .sort_values("n_edges", ascending=False)
           .reset_index())
    lad.to_csv(OUT / "delta_mz_ladders.tsv", sep="\t", index=False)

    # Figure: component size distribution before/after the artifact filter.
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    sizes_all = pd.Series([len(m) for m in members.values()])
    sizes_kept = df["size_in_analysis_matrix"]
    bins = np.arange(1, max(sizes_all.max(), 2) + 2) - 0.5
    ax.hist(sizes_all, bins=bins, color="#BBBBBB", label=f"network ({len(sizes_all)} components)")
    ax.hist(sizes_kept, bins=bins, color="#0072B2", alpha=0.85,
            label=f"in analysis matrix ({len(sizes_kept)} components)")
    ax.set_yscale("log")
    ax.set_xlabel("component size (features)")
    ax.set_ylabel("number of components (log)")
    ax.set_title("Molecular-family sizes shrink once isotope/adduct/ISF rows are removed", fontsize=10)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "component_sizes.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "component_sizes.pdf", bbox_inches="tight")
    plt.close(fig)

    # Markdown summary
    lines = ["# Molecular-network components", "",
             f"- network: **{len(edges)} edges, {n_nodes_all} nodes, {n_comp_all} components** "
             f"(largest {sizes_all.max()})",
             f"- after the artifact filter and analysis-matrix restriction: "
             f"**{len(df)} components**, largest **{int(sizes_kept.max())}**, "
             f"{int((sizes_kept >= 3).sum())} with >= 3 members",
             f"- background blank-clearing rate: " +
             ", ".join(f"{k} {v:.2%}" for k, v in bg.items()),
             ""]
    for sp_short in blank_ok:
        sig = df[(df[f"q_blank_enriched_{sp_short}"] < 0.05) & (df["size_in_analysis_matrix"] >= 3)]
        allb = df[df[f"all_blank_clearing_{sp_short}"]]
        lines += [f"## {sp_short}", "",
                  f"- components significantly blank-enriched (q<0.05, >=3 members): **{len(sig)}**",
                  f"- components ENTIRELY blank-clearing (>=3 members): **{len(allb)}**", ""]
        if len(sig):
            top = sig.head(10)[["component", "size_in_analysis_matrix",
                                f"n_blank_clearing_{sp_short}", f"q_blank_enriched_{sp_short}",
                                "dominant_npc_pathway", "example_structures"]]
            lines += ["| component | size | blank-clearing | q | dominant class | examples |",
                      "|---|---|---|---|---|---|"]
            for _, r in top.iterrows():
                lines.append(
                    f"| {int(r['component'])} | {int(r['size_in_analysis_matrix'])} | "
                    f"{int(r[f'n_blank_clearing_{sp_short}'])} | "
                    f"{r[f'q_blank_enriched_{sp_short}']:.2e} | "
                    f"{r['dominant_npc_pathway'] or '—'} | {r['example_structures'][:70] or '—'} |")
            lines.append("")
    lines += ["## DeltaMZ homologous steps (both ends surviving the artifact filter)", "",
              "| step | edges | median cosine |", "|---|---|---|"]
    for _, r in lad.iterrows():
        lines.append(f"| {r['delta_step']} | {int(r['n_edges'])} | "
                     f"{r['median_cosine']:.3f} |")
    lines += ["", f"Total classified edges: {int(lad['n_edges'].sum())} of {len(e)} "
                  f"intra-matrix edges ({lad['n_edges'].sum()/max(len(e),1):.1%}).", ""]
    (OUT / "component_summary.md").write_text("\n".join(lines))

    print(f"wrote {len(df)} components and {len(lad)} ladder classes to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
