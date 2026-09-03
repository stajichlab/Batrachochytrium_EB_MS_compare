#!/usr/bin/env python3
"""
Primary differential feature contrasts for the Bd/Bsal metabolomics table,
built on the collapsed life-stage vocabulary established by
analysis/ordination/scripts/build_ordination_table.py:

    stage_group = Zoospore | Developed (Sporangium + Mature collapsed)

The 30-way per-condition pairwise scan in analysis/differential_features/ is
the exploratory tier; this is the hypothesis tier with two contrast families
(stratified by species throughout -- media is a species confounder, Bd 1%
Tryptone vs Bsal 50% TGHL, so species are never merged):

  1. Life-stage (stage_group): within each species x matrix, Zoospore vs
     Developed -- the "zoospore vs the rest of the stages" question. Mature
     and Sporangium were near-identical (0 significant features in the
     within-matrix scan), so collapsing them buys power: liq n=10 vs 20,
     spore n=5 vs 10 per species.
        dendrobatidis {liq,spore} Zoospore_vs_Developed
        salamandrivorans {liq,spore} Zoospore_vs_Developed      (4 contrasts)

  2. Secreted-vs-cellular (matrix): within each species x stage_group, liq
     vs spore -- isolates features of the liquid-culture supernatant
     (extracellular/secreted fraction) vs the spore-pellet (cell-associated).
        dendrobatidis {Zoospore,Developed} liq_vs_spore
        salamandrivorans {Zoospore,Developed} liq_vs_spore      (4 contrasts)

Method per contrast (reuses the analysis/differential_features machinery --
Mann-Whitney U on prevalence-filtered, TSS-normalized peak areas; BH-FDR;
pseudocount-scaled log2 median fold-change).

Per-contrast outputs (analysis/differential_features_primary/<name>/):
  - differential_features.csv.gz : full ranked table.
  - volcano.png/.pdf, top_features.png/.pdf, top_features.tsv : as before.

Annotated / summary outputs (analysis/differential_features_primary/):
  - significant_annotated.tsv : every feature with q < --fdr across all 8
    primary contrasts, joined to analysis/sirius_annotation/
    sirius_annotations.tsv (NPC pathway/class, ClassyFire class, structure
    info) and flagged:
        liq_over_spore_log2fc : median(log2FC, liq vs spore) across that
          species' samples (stage-confounded context) -- positive = feature
          is more abundant in the liquid/supernatant among its parent stages,
          a hint it is enriched in the extracellular/secreted fraction.
        is_liq_enriched : liq_over_spore_log2fc >= +1 (the raw liq/spore
          direction, with NO media correction -- this was the pre-2026-09-02
          definition of is_secreted_candidate).
        passes_media_blank : feature is >= 2x its C_liq media-blank companion
          in at least one life stage for that species (background_subtraction.
          fungal_over_blank_ratio, the same filter the genome-bioactivity
          linkage pipeline applies at its Stage 5).
        is_secreted_candidate : is_liq_enriched AND passes_media_blank.
        bioactive : curated keyword hit suggesting high bioactivity.

  - significant_bioactive.tsv : subset of the above, bioactivity flagged;
    includes the per-species life-stage direction (up = higher in Zoospore).
  - primary_comparison_summary.tsv : one row per contrast, n per side,
    n_tested, n_significant, n_significant_annotated.

CAVEAT that motivates the media-blank term (added 2026-09-02): both growth
media are peptide-rich broths (Bd 1% tryptone = casein digest, Bsal 50%
TGHL). Media peptides are abundant in the `liq` supernatant and absent from
the washed `spore` pellet, so a raw liq-vs-spore contrast scores them as
maximally "secreted" -- the single largest false-positive class for this
goal. Requiring the feature to also exceed its own media blank removes them.
Without this term only ~9% (Bd) / ~20% (Bsal) of liq-enriched features are
distinguishable from medium.

Bioactivity keyword list is a curated heuristic (specialized-metabolism and
activity terms), applied to structure name + NPC pathway/class + ClassyFire
class; treat as a filter for manual curation, not a claim.

Usage:
    python3 scripts/differential_features_primary.py [--top-n 20] [--fdr 0.05] [--prevalence-min 0.10]
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

REPO = Path(__file__).resolve().parents[3]
LINKED = REPO / "analysis" / "ordination" / "linked_data"
OUT_ROOT = REPO / "analysis" / "differential_features_primary"
SIRIUS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"

# background_subtraction (+ its `paths` import) lives with the genome-bioactivity
# linkage pipeline; reuse it rather than reimplementing the blank contrast, so
# both pipelines apply an identical media-blank definition.
sys.path.insert(0, str(REPO / "analysis" / "genome_bioactivity_linkage" / "scripts"))
from background_subtraction import (  # noqa: E402
    fungal_over_blank_ratio,
    load_feature_intensities,
    load_metadata,
)

# The three sampled life stages; a feature only has to clear its blank in ONE
# of them to count (secretion is expected to be stage-specific).
BLANK_STAGES = ["Zoospore", "Sporangium", "Mature"]
BLANK_MIN_FC = 2.0


def media_blank_pass_ids(species_full: str) -> set[int]:
    """row_ids >= BLANK_MIN_FC x their C_liq companion blank in any life stage."""
    features = load_feature_intensities()
    meta = load_metadata()
    passing: set[int] = set()
    for stage in BLANK_STAGES:
        ratio = fungal_over_blank_ratio(
            features, meta, species=species_full, life_stage=stage, min_fc=BLANK_MIN_FC
        )
        passing |= set(ratio.loc[ratio["passes_background_filter"], "row_id"])
    print(
        f"[media-blank] {species_full}: {len(passing)} of {len(features)} features "
        f"clear {BLANK_MIN_FC:g}x their C_liq blank in >=1 life stage",
        file=sys.stderr,
    )
    return {int(r) for r in passing}

UP_COLOR = "#D55E00"  # vermillion -- higher in group A
DOWN_COLOR = "#0072B2"  # blue -- higher in group B
NS_COLOR = "#BBBBBB"

# Curated heuristic keyword set for "high bioactivity" flagging.
BIOACTIVE_RE = re.compile(
    r"(antibiotic|antimicrobial|antibacterial|antifungal|mycotoxin|phytoalexin|"
    r"phytotoxin|siderophore|quorum|pheromone|hormone|surfactant|toxin|"
    r"antineoplastic|anticancer|apoptosis|cytotoxic|mitogenic|immunomodulat|"
    r"ergosterol|sphingo|ceramide|alkaloid|terpenoid|polyketide|lipid)",
    re.IGNORECASE,
)

SPECIES_ORDER = ["Batrachochytrium dendrobatidis", "Batrachochytrium salamandrivorans"]
STAGE_PAIRS = [("Zoospore", "Developed")]


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.clip(ranked, 0, 1)
    return qvals


def tss_normalize(feat: pd.DataFrame, sample_ids: list[str], prevalence_min: float):
    mat = feat[sample_ids].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= prevalence_min
    kept_annot = feat.loc[keep, ["row_id", "mz", "rt"]].reset_index(drop=True)
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_ids, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = (mat / col_sums).T
    return kept_annot, mat_norm


def test_features(mat_a: np.ndarray, mat_b: np.ndarray) -> pd.DataFrame:
    n_features = mat_a.shape[1]
    pvals = np.empty(n_features)
    log2fc = np.empty(n_features)
    median_a = np.empty(n_features)
    median_b = np.empty(n_features)
    ustat = np.empty(n_features)

    for i in range(n_features):
        a, b = mat_a[:, i], mat_b[:, i]
        median_a[i], median_b[i] = np.median(a), np.median(b)
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        ustat[i], pvals[i] = u, p
        pseudocount = min(x[x > 0].min() if (x > 0).any() else 1e-12 for x in (a, b)) / 2
        log2fc[i] = np.log2((median_a[i] + pseudocount) / (median_b[i] + pseudocount))

    qvals = bh_fdr(pvals)
    return pd.DataFrame({
        "median_a": median_a, "median_b": median_b, "log2FC_a_over_b": log2fc,
        "U_stat": ustat, "p_value": pvals, "q_value": qvals,
    })


def savefig_multi(fig, png_path: Path):
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight")


def plot_volcano(df: pd.DataFrame, label_a: str, label_b: str, fdr: float, out_path: Path):
    sig_up = (df["q_value"] < fdr) & (df["log2FC_a_over_b"] > 0)
    sig_down = (df["q_value"] < fdr) & (df["log2FC_a_over_b"] < 0)
    ns = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(
        df.loc[ns, "log2FC_a_over_b"], -np.log10(df.loc[ns, "p_value"]),
        s=8, color=NS_COLOR, edgecolor="none", alpha=0.5, label=f"not significant (n={ns.sum()})",
    )
    ax.scatter(
        df.loc[sig_up, "log2FC_a_over_b"], -np.log10(df.loc[sig_up, "p_value"]),
        s=14, color=UP_COLOR, edgecolor="white", linewidth=0.3,
        label=f"higher in {label_a} (n={sig_up.sum()})",
    )
    ax.scatter(
        df.loc[sig_down, "log2FC_a_over_b"], -np.log10(df.loc[sig_down, "p_value"]),
        s=14, color=DOWN_COLOR, edgecolor="white", linewidth=0.3,
        label=f"higher in {label_b} (n={sig_down.sum()})",
    )
    ax.set_xlabel(f"log2 fold-change (median {label_a} / median {label_b})")
    ax.set_ylabel("-log10(p-value)")
    ax.set_title(f"Primary differential features: {label_a} vs {label_b}\n(FDR < {fdr:.0%} colored)")
    ax.axvline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_top_features(df: pd.DataFrame, mat_a, mat_b, label_a, label_b, out_path: Path, top_n: int):
    top = df.head(top_n)
    ncols = 4
    nrows = -(-top_n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows))
    axes = np.atleast_1d(axes).ravel()

    for ax, (_, row) in zip(axes, top.iterrows()):
        i = int(row.name)
        rng = np.random.default_rng(0)
        xa = 0 + rng.uniform(-0.12, 0.12, size=mat_a.shape[0])
        xb = 1 + rng.uniform(-0.12, 0.12, size=mat_b.shape[0])
        ax.scatter(xa, mat_a[:, i], s=14, color=UP_COLOR, edgecolor="none", alpha=0.8)
        ax.scatter(xb, mat_b[:, i], s=14, color=DOWN_COLOR, edgecolor="none", alpha=0.6)
        ax.set_xticks([0, 1])
        ax.set_xticklabels([label_a, label_b], fontsize=7)
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.set_title(
            f"row {row['row_id']:.0f}  m/z {row['mz']:.2f}\nq={row['q_value']:.1e}  log2FC={row['log2FC_a_over_b']:.2f}",
            fontsize=7,
        )
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(top):]:
        ax.set_visible(False)

    fig.suptitle(f"Top {top_n} primary differential features: {label_a} vs {label_b} (TSS proportion, symlog y)", fontsize=10)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def run_one(name: str, label_a: str, label_b: str, ids_a: list[str], ids_b: list[str],
            feat: pd.DataFrame, prevalence_min: float, fdr: float, top_n: int):
    kept_annot, mat_all = tss_normalize(feat, ids_a + ids_b, prevalence_min)
    mat_a, mat_b = mat_all[:len(ids_a)], mat_all[len(ids_a):]

    stats = test_features(mat_a, mat_b)
    result = pd.concat([kept_annot, stats], axis=1)
    result = result.assign(_abs=result["log2FC_a_over_b"].abs()).sort_values(
        ["q_value", "_abs"], ascending=[True, False]
    ).drop(columns="_abs")

    out_dir = OUT_ROOT / name
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "differential_features.csv.gz", index=False)
    plot_volcano(stats, label_a, label_b, fdr, out_dir / "volcano.png")
    plot_top_features(result, mat_a, mat_b, label_a, label_b, out_dir / "top_features.png", top_n)
    result.head(top_n).to_csv(out_dir / "top_features.tsv", sep="\t", index=False)

    n_sig = int((stats["q_value"] < fdr).sum())
    print(
        f"[{name}] {label_a} (n={len(ids_a)}) vs {label_b} (n={len(ids_b)}): "
        f"{len(stats)} features tested (prevalence >= {prevalence_min:.0%}), "
        f"{n_sig} significant at FDR < {fdr:.0%}",
        file=sys.stderr,
    )
    return out_dir.name, len(stats), n_sig


def flag_bioactive(row: pd.Series) -> bool:
    haystack = " | ".join(str(row.get(c, "")) for c in (
        "sirius_structure_name", "sirius_npc_pathway", "sirius_npc_class", "sirius_classyfire_class",
    ))
    return bool(haystack and BIOACTIVE_RE.search(haystack))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--prevalence-min", type=float, default=0.10)
    args = ap.parse_args()

    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz")
    sirius = pd.read_csv(SIRIUS, sep="\t")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    # Precompute per-feature stage-confounded liq-vs-spore direction, per species.
    liq_dir = {}
    for species in SPECIES_ORDER:
        sp_meta = meta[meta["species"] == species]
        ids_liq = sp_meta.loc[sp_meta["matrix"] == "liq", "sample_id"].tolist()
        ids_spore = sp_meta.loc[sp_meta["matrix"] == "spore", "sample_id"].tolist()
        annot, mat = tss_normalize(feat, ids_liq + ids_spore, args.prevalence_min)
        ml, ms = mat[:len(ids_liq)], mat[len(ids_liq):]
        fcs = np.empty(mat.shape[1])
        for i in range(mat.shape[1]):
            a, b = ml[:, i], ms[:, i]
            pc = min(x[x > 0].min() if (x > 0).any() else 1e-12 for x in (a, b)) / 2
            fcs[i] = np.log2((np.median(a) + pc) / (np.median(b) + pc))
        liq_dir[species] = pd.DataFrame({"row_id": annot["row_id"], "liq_over_spore_log2fc": fcs}).set_index("row_id")

    # Media-blank pass sets, keyed by the SHORT species name used in `combined`.
    blank_pass = {s.split()[-1]: media_blank_pass_ids(s) for s in SPECIES_ORDER}

    summary_rows = []
    sig_rows = []
    for species in SPECIES_ORDER:
        short = species.split()[-1]
        sp_meta = meta[meta["species"] == species]
        for matrix in ["liq", "spore"]:
            for stage_a, stage_b in STAGE_PAIRS:
                cond_a = f"{matrix}_{stage_a}"; cond_b = f"{matrix}_{stage_b}"
                ids_a = sp_meta.loc[sp_meta["condition_group"] == cond_a, "sample_id"].tolist()
                ids_b = sp_meta.loc[sp_meta["condition_group"] == cond_b, "sample_id"].tolist()
                name = f"{short}_{cond_a}_vs_{cond_b}"
                _, n_tested, n_sig = run_one(
                    name, cond_a, cond_b, ids_a, ids_b, feat,
                    args.prevalence_min, args.fdr, args.top_n,
                )
                summary_rows.append({
                    "comparison": name, "species": short, "family": "life_stage",
                    "group_a": cond_a, "group_b": cond_b, "n_a": len(ids_a), "n_b": len(ids_b),
                    "n_tested": n_tested, "n_significant": n_sig,
                })
                sub = pd.read_csv(OUT_ROOT / name / "differential_features.csv.gz")
                sub_sig = sub[sub["q_value"] < args.fdr].copy()
                if sub_sig.empty:
                    continue
                sub_sig["comparison"] = name
                sub_sig["species"] = short
                sub_sig["family"] = "life_stage"
                sub_sig["group_a"] = cond_a
                sub_sig["group_b"] = cond_b
                sub_sig["direction"] = np.where(
                    sub_sig["log2FC_a_over_b"] > 0, f"up_in_{cond_a}", f"up_in_{cond_b}"
                )
                # Secreted-candidate hint: liq/spore direction across this species.
                sub_sig = sub_sig.merge(liq_dir[species], on="row_id", how="left")
                sig_rows.append(sub_sig)

        for stage_group in ["Zoospore", "Developed"]:
            cond_a = f"liq_{stage_group}"; cond_b = f"spore_{stage_group}"
            ids_a = sp_meta.loc[sp_meta["condition_group"] == cond_a, "sample_id"].tolist()
            ids_b = sp_meta.loc[sp_meta["condition_group"] == cond_b, "sample_id"].tolist()
            name = f"{short}_{cond_a}_vs_{cond_b}"
            _, n_tested, n_sig = run_one(
                name, cond_a, cond_b, ids_a, ids_b, feat,
                args.prevalence_min, args.fdr, args.top_n,
            )
            summary_rows.append({
                "comparison": name, "species": short, "family": "secreted_vs_cellular",
                "group_a": cond_a, "group_b": cond_b, "n_a": len(ids_a), "n_b": len(ids_b),
                "n_tested": n_tested, "n_significant": n_sig,
            })
            sub = pd.read_csv(OUT_ROOT / name / "differential_features.csv.gz")
            sub_sig = sub[sub["q_value"] < args.fdr].copy()
            if sub_sig.empty:
                continue
            sub_sig["comparison"] = name
            sub_sig["species"] = short
            sub_sig["family"] = "secreted_vs_cellular"
            sub_sig["group_a"] = cond_a
            sub_sig["group_b"] = cond_b
            sub_sig["direction"] = np.where(
                sub_sig["log2FC_a_over_b"] > 0, f"up_in_{cond_a}", f"up_in_{cond_b}"
            )
            # Secreted-candidate hint: liq/spore direction across this species.
            # (Same merge as the life_stage family above -- without it,
            # liq_over_spore_log2fc/is_secreted_candidate are silently NaN/False
            # for every secreted_vs_cellular row, defeating the flag's purpose.)
            sub_sig = sub_sig.merge(liq_dir[species], on="row_id", how="left")
            sig_rows.append(sub_sig)

    if sig_rows:
        combined = pd.concat(sig_rows, ignore_index=True)
        combined["row_id"] = combined["row_id"].astype(int)
        sirius["row ID"] = sirius["row ID"].astype(int)
        combined = combined.merge(sirius, left_on="row_id", right_on="row ID", how="left", validate="many_to_one")
        # consolidate row_id float <-> sirius int key
        # Raw liq/spore direction (pre-2026-09-02 is_secreted_candidate).
        combined["is_liq_enriched"] = combined["liq_over_spore_log2fc"] >= 1.0
        # Media-blank term: the feature must also exceed its own C_liq blank.
        # Keyed on the SHORT species name, which is what `combined.species` holds.
        combined["passes_media_blank"] = [
            rid in blank_pass[sp] for rid, sp in zip(combined["row_id"], combined["species"])
        ]
        combined["is_secreted_candidate"] = (
            combined["is_liq_enriched"] & combined["passes_media_blank"]
        )
        bio = combined.apply(flag_bioactive, axis=1)
        combined["bioactive"] = bio
        combined = combined.sort_values(["species", "q_value"])
        combined.to_csv(OUT_ROOT / "significant_annotated.tsv", sep="\t", index=False)
        combined[combined["bioactive"]].to_csv(OUT_ROOT / "significant_bioactive.tsv", sep="\t", index=False)
        n_liq = int(combined["is_liq_enriched"].sum())
        n_sec = int(combined["is_secreted_candidate"].sum())
        print(
            f"wrote {len(combined)} significant (FDR < {args.fdr:.0%}) feature-rows "
            f"({int(bio.sum())} bioactivity-flagged) across {len(combined['comparison'].unique())} primary contrasts",
            file=sys.stderr,
        )
        print(
            f"  liq-enriched rows: {n_liq}; also clearing the media blank "
            f"(is_secreted_candidate): {n_sec} "
            f"({n_sec / n_liq:.1%} of liq-enriched)" if n_liq else "  liq-enriched rows: 0",
            file=sys.stderr,
        )

    summary = pd.DataFrame(summary_rows).sort_values(["family", "species", "n_significant"], ascending=[True, True, False])
    summary.to_csv(OUT_ROOT / "primary_comparison_summary.tsv", sep="\t", index=False)
    print(f"wrote {len(summary)} primary comparisons to {OUT_ROOT}", file=sys.stderr)


if __name__ == "__main__":
    main()
