#!/usr/bin/env python3
"""
Pairwise differential MS2 feature testing between the 6 (matrix,
life_stage) states, within each species -- "what actually differs" behind
the ordination in analysis/ordination/ (F-001: matrix dominates; this asks
which features drive that and every other pairwise contrast).

Reuses analysis/ordination/linked_data/{sample_metadata.csv,
feature_abundance.csv.gz} directly (same 90-sample table the PCoA scripts
use) rather than re-deriving it.

Method, per feature, per comparison (mirrors the Rhodotorula sibling
project's scripts/differential_features_by_species.py):
  1. Prevalence filter (>=10% of samples in the two groups being compared)
     + total-sum-scaling (TSS) -- no power transform here (that's for
     stabilizing a Bray-Curtis distance matrix; fold-change is reported on
     the directly-interpretable TSS proportions instead).
  2. Mann-Whitney U test (rank-based; robust to n=5-10/group and to the
     heavy right-skew typical of peak-area data).
  3. Effect size = log2 fold-change of group medians, with a pseudocount
     that scales with the feature (half the smallest nonzero TSS
     proportion observed for that feature across both groups).
  4. Benjamini-Hochberg FDR across all tested features (no `statsmodels`
     installed in this repo's pixi env -- implemented directly, same
     order-statistic formula `statsmodels.stats.multitest.multipletests`
     uses).

Comparisons: every pairwise contrast among the 6 (matrix, life_stage)
states, within each species (15 pairs x 2 species = 30 total) --
matches analysis/ordination's condition vocabulary.

Output per comparison (analysis/differential_features/<species>_<condA>_vs_<condB>/):
  - differential_features.csv.gz : full ranked table, one row per tested
    feature (row_id, mz, rt, group medians, log2FC, U-stat, p-value,
    q-value), sorted by q-value then |log2FC|.
  - volcano.png/.pdf   : log2FC vs -log10(p), FDR-significant features
    colored by direction (dataviz convention: null test is the FDR
    threshold itself here, not a permutation test -- see caveats in
    analysis/differential_features/README.md).
  - top_features.png/.pdf : per-sample TSS abundance (small multiples) for
    the top TOP_N features by q-value.
  - top_features.tsv  : the same top-TOP_N features as a plain table.

Usage:
    python3 scripts/differential_features.py [--top-n 20] [--fdr 0.05] [--prevalence-min 0.10]
"""
import argparse
import itertools
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
OUT_ROOT = REPO / "analysis" / "differential_features"

CONDITION_ORDER = [
    "liq_Zoospore", "liq_Sporangium", "liq_Mature",
    "spore_Zoospore", "spore_Sporangium", "spore_Mature",
]
UP_COLOR = "#D55E00"  # vermillion -- higher in group A
DOWN_COLOR = "#0072B2"  # blue -- higher in group B
NS_COLOR = "#BBBBBB"


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR, same order-statistic formula as
    statsmodels.stats.multitest.multipletests(method='fdr_bh')."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    qvals = np.empty(n)
    qvals[order] = np.clip(ranked, 0, 1)
    return qvals


def tss_normalize(feat: pd.DataFrame, sample_ids: list[str], prevalence_min: float):
    """Prevalence-filter then TSS-normalize. Returns (kept feature-id
    DataFrame with a fresh index, samples x features ndarray)."""
    mat = feat[sample_ids].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= prevalence_min
    kept_annot = feat.loc[keep, ["row_id", "mz", "rt"]].reset_index(drop=True)
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_ids, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = (mat / col_sums).T  # samples x features
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
    ax.set_title(f"Differential MS2 features: {label_a} vs {label_b}\n(FDR < {fdr:.0%} colored)")
    ax.axvline(0, color="#DDDDDD", linewidth=0.8, zorder=0)
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_top_features(df: pd.DataFrame, mat_a, mat_b, label_a, label_b, out_path: Path, top_n: int):
    """df must be sorted by significance and have a fresh 0..n-1 RangeIndex
    matching mat_a/mat_b's columns."""
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

    fig.suptitle(f"Top {top_n} differential features: {label_a} vs {label_b} (TSS proportion, symlog y)", fontsize=10)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def run_one(species_short: str, cond_a: str, cond_b: str, meta: pd.DataFrame, feat: pd.DataFrame,
            prevalence_min: float, fdr: float, top_n: int):
    ids_a = meta.loc[meta["condition"] == cond_a, "sample_id"].tolist()
    ids_b = meta.loc[meta["condition"] == cond_b, "sample_id"].tolist()
    kept_annot, mat_all = tss_normalize(feat, ids_a + ids_b, prevalence_min)
    mat_a, mat_b = mat_all[:len(ids_a)], mat_all[len(ids_a):]

    stats = test_features(mat_a, mat_b)
    result = pd.concat([kept_annot, stats], axis=1)
    result = result.assign(_abs=result["log2FC_a_over_b"].abs()).sort_values(
        ["q_value", "_abs"], ascending=[True, False]
    ).drop(columns="_abs")

    out_dir = OUT_ROOT / f"{species_short}_{cond_a}_vs_{cond_b}"
    out_dir.mkdir(parents=True, exist_ok=True)
    result.to_csv(out_dir / "differential_features.csv.gz", index=False)
    plot_volcano(stats, cond_a, cond_b, fdr, out_dir / "volcano.png")
    plot_top_features(result, mat_a, mat_b, cond_a, cond_b, out_dir / "top_features.png", top_n)
    result.head(top_n).to_csv(out_dir / "top_features.tsv", sep="\t", index=False)

    n_sig = int((stats["q_value"] < fdr).sum())
    print(
        f"[{species_short}] {cond_a} (n={len(ids_a)}) vs {cond_b} (n={len(ids_b)}): "
        f"{len(stats)}/{len(feat)} features tested (prevalence >= {prevalence_min:.0%}), "
        f"{n_sig} significant at FDR < {fdr:.0%}",
        file=sys.stderr,
    )
    return out_dir.name, len(stats), n_sig


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-n", type=int, default=20)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--prevalence-min", type=float, default=0.10)
    args = ap.parse_args()

    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz")
    OUT_ROOT.mkdir(parents=True, exist_ok=True)

    summary_rows = []
    for species in sorted(meta["species"].unique()):
        short = species.split()[-1]
        sp_meta = meta[meta["species"] == species]
        conds_present = [c for c in CONDITION_ORDER if (sp_meta["condition"] == c).any()]
        for cond_a, cond_b in itertools.combinations(conds_present, 2):
            name, n_tested, n_sig = run_one(
                short, cond_a, cond_b, sp_meta, feat, args.prevalence_min, args.fdr, args.top_n
            )
            summary_rows.append({
                "comparison": name, "species": short, "condition_a": cond_a, "condition_b": cond_b,
                "n_tested": n_tested, "n_significant": n_sig,
            })

    summary = pd.DataFrame(summary_rows).sort_values("n_significant", ascending=False)
    summary.to_csv(OUT_ROOT / "comparison_summary.csv", index=False)
    print(f"wrote {len(summary)} comparisons to {OUT_ROOT}", file=sys.stderr)


if __name__ == "__main__":
    main()
