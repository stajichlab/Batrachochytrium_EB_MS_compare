#!/usr/bin/env python3
"""
Bray-Curtis PCoA ordinations of the Bd/Bsal MS2 feature-abundance table
(90 analysis samples: 2 species x 2 matrix (liq/spore) x 3 life_stage
(Zoospore/Sporangium/Mature)), following the dataviz convention pack
(.living/conventions/dataviz/analysis-conventions.md): variance-explained
axis labels, colorblind-safe (Okabe-Ito) palette, vector output alongside
PNG.

Preprocessing (independent per sample subset, since prevalence depends on
which samples are included -- see Rhodotorula sibling-project pattern in
scripts/pcoa_ms_features.py):
  1. Prevalence filter: keep features with peak area > 0 in >= 10% of the
     samples in that subset.
  2. Total-sum-scaling (TSS) per sample, then a fourth-root power transform
     (variance-stabilizing, stays non-negative so Bray-Curtis is still
     defined -- a log transform would produce negatives).
  3. Bray-Curtis distance + classical (Torgerson/Gower) PCoA; variance
     explained is of the *positive* eigenvalues only (Bray-Curtis is
     non-Euclidean so classical PCoA has small negative eigenvalues by
     construction; the convention here matches the Rhodotorula reference
     script and is noted in each figure title).

Two tiers of figures:
  1. `pcoa_all.png`     -- every analysis sample, color=species,
                            shape=matrix (liq/spore).
  2. `by_species/`      -- one ordination per species, color=life_stage,
                            shape=matrix.

Per-pair state contrasts (matrix/life_stage pairwise) are handled instead
by the feature-level differential/volcano pipeline in
analysis/differential_features/ -- a PCoA restricted to just two groups at
a time is not an especially informative view (see 2026-08-18 session note),
so that tier was dropped from here.

Usage:
    python3 scripts/pcoa_ordination.py
"""
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

REPO = Path(__file__).resolve().parents[3]
LINKED = REPO / "analysis" / "ordination" / "linked_data"
FIG_DIR = REPO / "analysis" / "ordination" / "figures"
PREVALENCE_MIN = 0.10

# Okabe-Ito colorblind-safe categorical palette.
SPECIES_COLORS = {
    "Batrachochytrium dendrobatidis": "#0072B2",
    "Batrachochytrium salamandrivorans": "#D55E00",
}
LIFESTAGE_COLORS = {
    "Zoospore": "#0072B2",
    "Sporangium": "#E69F00",
    "Mature": "#009E73",
}
MATRIX_MARKERS = {"liq": "o", "spore": "^"}


def prep_matrix(feat: pd.DataFrame, sample_ids: list[str]):
    mat = feat[sample_ids].to_numpy(dtype=float)
    prevalence = (mat > 0).mean(axis=1)
    keep = prevalence >= PREVALENCE_MIN
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, ok in zip(sample_ids, col_sums == 0) if ok]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    mat_norm = (mat / col_sums) ** 0.25
    return mat_norm.T, int(keep.sum()), len(keep)


def classical_pcoa(dist: np.ndarray):
    n = dist.shape[0]
    d2 = dist**2
    centering = np.eye(n) - np.ones((n, n)) / n
    b = -0.5 * centering @ d2 @ centering
    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]
    pos = eigvals > 1e-8 * eigvals[0]
    coords = eigvecs[:, pos] * np.sqrt(eigvals[pos])
    prop = eigvals[pos] / eigvals[pos].sum()
    return coords, prop, eigvals


def run_ordination(name: str, sample_ids: list[str], feat: pd.DataFrame):
    mat, n_kept, n_total = prep_matrix(feat, sample_ids)
    dist = squareform(pdist(mat, metric="braycurtis"))
    coords, prop, eigvals = classical_pcoa(dist)
    n_neg = int((eigvals < -1e-8).sum())
    print(
        f"[{name}] {len(sample_ids)} samples, {n_kept}/{n_total} features kept "
        f"(prevalence >= {PREVALENCE_MIN:.0%}), axis1={prop[0]:.1%} axis2={prop[1]:.1%} "
        f"of positive-eigenvalue variance ({n_neg} negative eigenvalues excluded)",
        file=sys.stderr,
    )
    axes = pd.DataFrame({"sample_id": sample_ids, "PCoA1": coords[:, 0], "PCoA2": coords[:, 1]})
    return axes, prop, dist


def savefig_multi(fig, png_path: Path):
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight")


def plot_all(axes, meta, prop, out_path):
    df = axes.merge(meta[["sample_id", "species", "matrix"]], on="sample_id")
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for species, color in SPECIES_COLORS.items():
        for matrix, marker in MATRIX_MARKERS.items():
            sub = df[(df["species"] == species) & (df["matrix"] == matrix)]
            if sub.empty:
                continue
            ax.scatter(
                sub["PCoA1"], sub["PCoA2"], marker=marker, s=32, color=color,
                edgecolor="white", linewidth=0.4,
                label=f"{species.split()[-1]} / {matrix} (n={len(sub)})",
            )
    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    ax.set_title("Bray-Curtis PCoA -- all analysis samples (n=90)")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def plot_by_species(axes, meta, prop, species, out_path):
    df = axes.merge(meta[["sample_id", "life_stage", "matrix"]], on="sample_id")
    fig, ax = plt.subplots(figsize=(7, 6))
    for stage, color in LIFESTAGE_COLORS.items():
        for matrix, marker in MATRIX_MARKERS.items():
            sub = df[(df["life_stage"] == stage) & (df["matrix"] == matrix)]
            if sub.empty:
                continue
            ax.scatter(
                sub["PCoA1"], sub["PCoA2"], marker=marker, s=36, color=color,
                edgecolor="white", linewidth=0.4,
                label=f"{stage} / {matrix} (n={len(sub)})",
            )
    ax.set_xlabel(f"PCoA1 ({prop[0]:.1%})")
    ax.set_ylabel(f"PCoA2 ({prop[1]:.1%})")
    ax.set_title(f"Bray-Curtis PCoA -- {species}")
    ax.legend(frameon=False, fontsize=8, loc="center left", bbox_to_anchor=(1.0, 0.5))
    ax.set_aspect("equal", adjustable="datalim")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    savefig_multi(fig, out_path)
    plt.close(fig)


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    (FIG_DIR / "by_species").mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz")

    all_ids = meta["sample_id"].tolist()
    axes_all, prop_all, _ = run_ordination("all", all_ids, feat)
    axes_all.to_csv(FIG_DIR / "pcoa_axes_all.csv", index=False)
    plot_all(axes_all, meta, prop_all, FIG_DIR / "pcoa_all.png")

    for species in sorted(meta["species"].unique()):
        ids = meta.loc[meta["species"] == species, "sample_id"].tolist()
        short = species.split()[-1]
        axes_sp, prop_sp, _ = run_ordination(short, ids, feat)
        axes_sp.to_csv(FIG_DIR / "by_species" / f"pcoa_axes_{short}.csv", index=False)
        plot_by_species(
            axes_sp, meta, prop_sp, species, FIG_DIR / "by_species" / f"pcoa_{short}.png"
        )

    print(f"wrote ordination figures to {FIG_DIR}", file=sys.stderr)


if __name__ == "__main__":
    main()
