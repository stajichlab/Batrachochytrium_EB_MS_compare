#!/usr/bin/env python3
"""Ordinal life-stage trend test across the THREE sampled stages.

Motivation (learnings L-21, 2026-09-02)
---------------------------------------
The primary tier (`differential_features_primary.py`) collapses
Sporangium+Mature into a single `Developed` stage. That collapse was
justified by the claim that "every within-matrix Sporangium-vs-Mature
contrast was 0-significant" -- which is FALSE for Bd's spore fraction:
`dendrobatidis_spore_Sporangium_vs_spore_Mature` = 5,507 significant of
21,816 tested in the 30-way scan, against 0 for Bsal spore.

Note the "Bsal spore is 2-state" reading of that 0 does NOT hold up. At
n=5v5 the minimum attainable two-sided Mann-Whitney p is 2/C(10,5)=0.0079,
so BH at q<0.05 needs ~16% of features at the floor simultaneously. Bsal's
Sporangium-vs-Mature separation is enriched 8.8x over the analytic null and
Bd's 26.6x -- both are real; only Bd's clears a discontinuous threshold.
Both species have a 3-state spore trajectory; the collapse loses signal in
both, which is the whole reason this trend tier exists.

Collapsing therefore discards real Bd signal. The 30-way pairwise scan in
`analysis/differential_features/` already covers every stage PAIR, so
re-running pairs here would duplicate it. What no tier currently provides is
a test of MONOTONIC progression across the ordered stages -- exactly the
structure a developmental trajectory has, and exactly what a two-group
collapse destroys.

Method
------
Per species x matrix (4 strata), on the same prevalence-filtered,
TSS-normalized matrix the other tiers use:

  stage rank: Zoospore=0, Sporangium=1, Mature=2 (ordered by timepoint
  8/48/96 h, so the rank is the real time axis, not an arbitrary coding)

  Spearman rho between each feature's normalized abundance and stage rank,
  two-sided, BH-FDR across features within the stratum.

Spearman (not Jonckheere-Terpstra) because with 3 x-levels the two are
equivalent up to a monotone transform; it is rank-based and gives a signed
effect size (rho: + = rises toward Mature, - = falls toward Zoospore).
Significance uses a LABEL-PERMUTATION null (see spearman_permutation), not
scipy's t-approximation, which is invalid under 3 tied x-levels.
n = 15 for spore strata (5 per stage) and n = 30 for liq strata (10 per
stage) -- ALL FUNGAL as of 2026-09-02; the liq strata previously contained
15 fungal + 15 sterile media blanks, which diluted every liq result.

Outputs (analysis/differential_features_primary/lifestage_trend/):
  <species>_<matrix>_trend.tsv  : every tested feature, rho/p/q, direction,
                                  SIRIUS identity, media-blank status
  trend_summary.tsv             : one row per stratum
  trend_<species>_<matrix>.png  : top monotonic features across the 3 stages

Caveats
-------
- n=15 (spore) bounds the achievable Spearman p-value; treat spore-stratum
  q-values as screening, not confirmation.
- These are the SAME samples the pairwise tiers use, so trend hits are not
  independent evidence from the pairwise hits -- this is a different
  question (monotonic progression) on shared data, not a replication.
- `liq_blank_status` is meaningful ONLY in the liq strata, where it tells a
  medium being consumed apart from fungal production. There is NO spore-pellet
  process blank in this design, so in the spore strata the column is written as
  "n/a (no spore blank exists)" rather than silently reusing the supernatant
  blank. Doing the latter asks whether a cell-associated compound exceeds a
  supernatant blank, and would have discarded the most unambiguously fungal
  features -- those never detected in any liq sample at all.

Usage:
    python3 scripts/lifestage_trend.py [--fdr 0.05] [--prevalence-min 0.10] [--top-n 12]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

REPO = Path(__file__).resolve().parents[3]
LINKED = REPO / "analysis" / "ordination" / "linked_data"
OUT = REPO / "analysis" / "differential_features_primary" / "lifestage_trend"
SIRIUS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"

sys.path.insert(0, str(REPO / "analysis" / "genome_bioactivity_linkage" / "scripts"))
from background_subtraction import (  # noqa: E402
    fungal_over_blank_ratio,
    load_feature_intensities,
    load_metadata,
)

SPECIES_ORDER = ["Batrachochytrium dendrobatidis", "Batrachochytrium salamandrivorans"]
# Ordered by timepoint (8 / 48 / 96 h) -- the rank IS the time axis.
STAGE_RANK = {"Zoospore": 0, "Sporangium": 1, "Mature": 2}
BLANK_STAGES = list(STAGE_RANK)
BLANK_MIN_FC = 2.0

UP_COLOR = "#D55E00"
DOWN_COLOR = "#0072B2"


def spearman_permutation(mat: np.ndarray, ranks: np.ndarray, n_perm: int, seed: int):
    """Spearman rho with a LABEL-PERMUTATION null instead of scipy's t-approximation.

    Why (corrected 2026-09-02): with only 3 distinct x-levels and heavy ties,
    scipy's asymptotic t-approximation is invalid and strongly
    anti-conservative -- for a perfectly ordered feature it returns 1.13e-7
    (n=15) and 6.09e-15 (n=30) against exact permutation floors of 2.64e-6 and
    3.60e-13, i.e. 23x and 59x too small. The first version of this script
    emitted ~1,500 p-values BELOW the exact attainable floor, which cannot
    exist. Discovery sets survived a 5-shuffle check, so the sets were real,
    but the q-values were not interpretable as FDR.

    Spearman itself is kept: with 3 x-levels it is equivalent to
    Jonckheere-Terpstra up to a monotone transform, so the test statistic was
    never the problem -- only its null distribution.

    Implementation: feature ranks are invariant under label permutation, so
    rho reduces to a Pearson correlation between standardized feature ranks
    and the standardized stage vector, and each permutation is a single
    matrix-vector product. The null is pooled across features (all features
    see the same shuffled labels), giving n_perm x n_features null values for
    empirical p resolution far finer than n_perm alone.
    """
    n_samp, n_feat = mat.shape
    # Rank within each feature (ties -> average), then standardize.
    fr = np.apply_along_axis(rankdata, 0, mat)
    fr -= fr.mean(axis=0)
    denom = np.sqrt((fr ** 2).sum(axis=0))
    constant = denom == 0
    denom[constant] = 1.0
    fr /= denom

    def rho_for(x: np.ndarray) -> np.ndarray:
        xs = x - x.mean()
        nrm = np.sqrt((xs ** 2).sum())
        return fr.T @ (xs / nrm)

    rho = rho_for(ranks)
    rho[constant] = 0.0

    rng = np.random.default_rng(seed)
    perm = ranks.copy()
    null = np.empty(n_perm * n_feat, dtype=np.float32)
    for k in range(n_perm):
        rng.shuffle(perm)
        null[k * n_feat:(k + 1) * n_feat] = np.abs(rho_for(perm))
    null.sort()

    # Empirical two-sided p: fraction of the pooled null at least as extreme.
    exceed = len(null) - np.searchsorted(null, np.abs(rho), side="left")
    pval = (exceed + 1) / (len(null) + 1)
    pval[constant] = 1.0
    return rho, pval, int(constant.sum())


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def tss_normalize(feat: pd.DataFrame, sample_ids: list[str], prevalence_min: float):
    mat = feat[sample_ids].to_numpy(dtype=float)
    keep = (mat > 0).mean(axis=1) >= prevalence_min
    annot = feat.loc[keep, ["row_id", "mz", "rt"]].reset_index(drop=True)
    mat = mat[keep]
    col_sums = mat.sum(axis=0)
    if (col_sums == 0).any():
        empty = [s for s, z in zip(sample_ids, col_sums == 0) if z]
        sys.exit(f"sample(s) with zero total abundance after filtering: {empty}")
    return annot, (mat / col_sums).T


def media_blank_pass_ids(species_full: str) -> set[int]:
    features = load_feature_intensities()
    meta = load_metadata()
    passing: set[int] = set()
    for stage in BLANK_STAGES:
        r = fungal_over_blank_ratio(
            features, meta, species=species_full, life_stage=stage, min_fc=BLANK_MIN_FC
        )
        passing |= set(r.loc[r["passes_background_filter"], "row_id"])
    return {int(x) for x in passing}


def plot_top(df, mat, ranks, species_short, matrix, out_path: Path, top_n: int):
    top = df.head(top_n)
    if top.empty:
        return
    ncols = 4
    nrows = -(-len(top) // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 2.6 * nrows))
    axes = np.atleast_1d(axes).ravel()
    rng = np.random.default_rng(0)
    labels = ["Zoo", "Spor", "Mat"]
    for ax, (_, row) in zip(axes, top.iterrows()):
        i = int(row["_col"])
        colour = UP_COLOR if row["rho"] > 0 else DOWN_COLOR
        ax.scatter(ranks + rng.uniform(-0.12, 0.12, size=len(ranks)), mat[:, i],
                   s=16, color=colour, edgecolor="none", alpha=0.8)
        med = [np.median(mat[ranks == r, i]) for r in (0, 1, 2)]
        ax.plot([0, 1, 2], med, color="#333333", linewidth=1.2, zorder=3)
        ax.set_xticks([0, 1, 2])
        ax.set_xticklabels(labels, fontsize=7)
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.set_title(f"row {int(row['row_id'])}  m/z {row['mz']:.2f}\n"
                     f"rho={row['rho']:+.2f}  q={row['q_value']:.1e}", fontsize=7)
        ax.spines[["top", "right"]].set_visible(False)
    for ax in axes[len(top):]:
        ax.set_visible(False)
    fig.suptitle(f"Monotonic life-stage trend: {species_short} {matrix} "
                 f"(TSS proportion, symlog y; line = per-stage median)", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    fig.savefig(out_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--fdr", type=float, default=0.05)
    ap.add_argument("--prevalence-min", type=float, default=0.10)
    ap.add_argument("--top-n", type=int, default=12)
    ap.add_argument("--n-perm", type=int, default=1000,
                    help="label permutations for the null (default 1000)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz")
    sirius = pd.read_csv(SIRIUS, sep="\t")
    sirius["row ID"] = sirius["row ID"].astype(int)

    summary = []
    for species in SPECIES_ORDER:
        short = species.split()[-1]
        blank_ok = media_blank_pass_ids(species)
        sp_meta = meta[meta["species"] == species]
        for matrix in ["liq", "spore"]:
            sub = sp_meta[sp_meta["matrix"] == matrix].copy()
            sub["_rank"] = sub["life_stage"].map(STAGE_RANK)
            if sub["_rank"].isna().any():
                sys.exit(f"unmapped life_stage in {short}/{matrix}: "
                         f"{sorted(sub.loc[sub['_rank'].isna(), 'life_stage'].unique())}")
            sub = sub.sort_values("_rank")
            ids = sub["sample_id"].tolist()
            ranks = sub["_rank"].to_numpy(dtype=float)
            if len(set(ranks)) != 3:
                sys.exit(f"{short}/{matrix} does not have all three stages")

            annot, mat = tss_normalize(feat, ids, args.prevalence_min)
            n_feat = mat.shape[1]
            rho, pval, n_const = spearman_permutation(
                mat, ranks, n_perm=args.n_perm, seed=args.seed
            )
            q = bh_fdr(pval)
            print(f"  ({n_const} features constant within stratum -> rho=0, p=1)",
                  file=sys.stderr)

            df = annot.copy()
            df["_col"] = np.arange(n_feat)
            df["rho"] = rho
            df["p_value"] = pval
            df["q_value"] = q
            df["direction"] = np.where(rho > 0, "rises_to_Mature", "falls_to_Mature")
            df["species"] = short
            df["matrix"] = matrix
            df["n_samples"] = len(ids)
            df["row_id"] = df["row_id"].astype(int)
            # Only meaningful for liq; no spore-pellet process blank exists.
            df["liq_blank_status"] = (
                df["row_id"].isin(blank_ok).map({True: "clears_blank", False: "at_or_below_blank"})
                if matrix == "liq" else "n/a (no spore blank exists)"
            )
            df = df.merge(sirius, left_on="row_id", right_on="row ID",
                          how="left", validate="many_to_one")
            df = df.sort_values(["q_value", "rho"], key=lambda s: s if s.name == "q_value" else -s.abs())

            n_sig = int((df["q_value"] < args.fdr).sum())
            n_sig_blank = (
                int(((df["q_value"] < args.fdr) & (df["liq_blank_status"] == "clears_blank")).sum())
                if matrix == "liq" else -1
            )
            df.drop(columns=["_col"]).to_csv(OUT / f"{short}_{matrix}_trend.tsv",
                                             sep="\t", index=False)
            plot_top(df[df["q_value"] < args.fdr], mat, ranks, short, matrix,
                     OUT / f"trend_{short}_{matrix}.png", args.top_n)

            summary.append({
                "species": short, "matrix": matrix, "n_samples": len(ids),
                "n_tested": n_feat, "n_significant": n_sig,
                "n_significant_blank_clearing": n_sig_blank,  # -1 = n/a (spore)
                "n_rises_to_Mature": int(((df["q_value"] < args.fdr) & (df["rho"] > 0)).sum()),
                "n_falls_to_Mature": int(((df["q_value"] < args.fdr) & (df["rho"] < 0)).sum()),
            })
            blank_note = (f"{n_sig_blank} also clear the media blank"
                          if matrix == "liq" else "no spore blank exists")
            print(f"[{short}/{matrix}] n={len(ids)}, {n_feat} tested, "
                  f"{n_sig} monotonic at FDR<{args.fdr:.0%} ({blank_note})",
                  file=sys.stderr)

    pd.DataFrame(summary).to_csv(OUT / "trend_summary.tsv", sep="\t", index=False)
    print(f"wrote {len(summary)} strata to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
