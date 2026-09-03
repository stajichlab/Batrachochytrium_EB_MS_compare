#!/usr/bin/env python3
"""Threshold-free effect measure for the small-n contrasts: complete-separation enrichment.

Why this exists
---------------
Every contrast in this project has n=5 or n=5-vs-10 per side. At those sizes
the Mann-Whitney U p-value has a hard floor:

    minimum attainable two-sided p = 2 / C(n1+n2, n1)

    n1=5,  n2=5  -> 2/252   = 7.94e-3
    n1=5,  n2=10 -> 2/3003  = 6.66e-4
    n1=10, n2=10 -> 2/184756 = 1.08e-5

A feature cannot be individually significant after BH correction unless a
LARGE NUMBER of features sit at that floor simultaneously: BH declares the
k-th smallest p significant when p_(k) <= q*k/m, so with every hit pinned at
p_min the requirement is

    k >= p_min * m / q

i.e. at n=5v5 and q=0.05, k >= 0.159*m -- about 16% of all tested features
must separate perfectly AT ONCE or nothing is called at all.

That makes "n_significant" a step function of the feature universe, not a
measure of biological effect. `dendrobatidis_spore_Sporangium_vs_spore_Mature`
illustrates it: 5,507 significant on the 38,547-feature table, 2,583 on the
artifact-filtered 25,157-feature table, and 0 in an intermediate run that used
scipy's asymptotic null (whose 5v5 floor of 1.219e-2 raises the required
simultaneous-separation count by 1.5x). The underlying separation count barely
moves across all three -- it sits at ~24x the null throughout. Only the
denominator and the p-floor move. Reporting any of those counts as "the
developmental signal" is an artifact of those two quantities.

What this script reports instead
--------------------------------
For each contrast, the number of features showing COMPLETE SEPARATION
(every value in one group strictly exceeds every value in the other), which
is exactly the event that attains p_min, compared against its analytic null
expectation:

    expected   = n_tested * 2 / C(n1+n2, n1)
    enrichment = observed / expected
    perm_p     = exact label-permutation p over every distinct relabelling

The null is a LABEL PERMUTATION, not a binomial. A binomial over features
would assume feature independence, which is exactly what the ordination
disproves, and it returned p = 0.0 for 10 of 12 contrasts. Permuting sample
labels keeps the between-feature correlation intact. Note the analytic
2/C(n,k) expectation is CONSERVATIVE here -- zero-inflation ties break
separation, so the permutation median runs below it and the enrichment ratios
are understated rather than inflated.

Enrichment is denominator-stable and comparable across contrasts, and it
answers "is there ordered signal here" without a discontinuous threshold.
Ties (shared zeros) break separation, so the measure is conservative on
sparse features.

This does NOT replace the per-feature FDR tables -- it is the honest summary
statistic for whether a contrast carries signal at all, and the correct basis
for statements like "Sporangium and Mature are indistinguishable".

Usage:
    python3 scripts/separation_enrichment.py
Outputs:
    analysis/differential_features/separation_enrichment.tsv
    analysis/differential_features/separation_enrichment.png/.pdf
"""
from __future__ import annotations

import sys
from math import comb
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from itertools import combinations

REPO = Path(__file__).resolve().parents[3]
LINKED = REPO / "analysis" / "ordination" / "linked_data"
OUT_DIR = REPO / "analysis" / "differential_features"

PREVALENCE_MIN = 0.10
STAGES = ["Zoospore", "Sporangium", "Mature"]
SPECIES = ["Batrachochytrium dendrobatidis", "Batrachochytrium salamandrivorans"]


def tss(feat: pd.DataFrame, ids: list[str]):
    mat = feat[ids].to_numpy(dtype=float)
    keep = (mat > 0).mean(axis=1) >= PREVALENCE_MIN
    mat = mat[keep]
    cs = mat.sum(axis=0)
    if (cs == 0).any():
        sys.exit(f"zero-total sample among {ids}")
    return (mat / cs).T


def n_complete_separations(a: np.ndarray, b: np.ndarray) -> int:
    """Features where min(A) > max(B) or max(A) < min(B) (strict, ties break it)."""
    return int(((a.min(axis=0) > b.max(axis=0)) | (a.max(axis=0) < b.min(axis=0))).sum())


def permutation_separations(mat: np.ndarray, n_a: int) -> np.ndarray:
    """Complete-separation counts under every distinct relabelling of the samples.

    At n=5 v 5 there are C(10,5)=252 assignments, but the separation count is
    identical for a split and its complement (separation is direction-agnostic
    here), so 126 distinct values -- the attainable p floor is 1/127.
    """
    n = mat.shape[0]
    seen, out = set(), []
    for idx in combinations(range(n), n_a):
        key = frozenset(idx)
        comp_key = frozenset(set(range(n)) - set(idx))
        if key in seen or comp_key in seen:
            continue
        seen.add(key)
        mask = np.zeros(n, dtype=bool)
        mask[list(idx)] = True
        out.append(n_complete_separations(mat[mask], mat[~mask]))
    return np.array(out)


def main() -> None:
    meta = pd.read_csv(LINKED / "sample_metadata.csv")
    feat = pd.read_csv(LINKED / "feature_abundance.csv.gz")
    if meta["is_C_companion"].any():
        sys.exit("media blanks present in sample_metadata.csv -- rebuild linked_data first")

    rows = []
    for species in SPECIES:
        short = species.split()[-1]
        sm = meta[meta["species"] == species]
        for matrix in ["liq", "spore"]:
            for i, sa in enumerate(STAGES):
                for sb in STAGES[i + 1:]:
                    ia = sm[(sm["matrix"] == matrix) & (sm["life_stage"] == sa)]["sample_id"].tolist()
                    ib = sm[(sm["matrix"] == matrix) & (sm["life_stage"] == sb)]["sample_id"].tolist()
                    mat = tss(feat, ia + ib)
                    A, B = mat[:len(ia)], mat[len(ia):]
                    m = mat.shape[1]
                    obs = n_complete_separations(A, B)
                    p_min = 2 / comb(len(ia) + len(ib), len(ia))
                    exp = m * p_min
                    # LABEL-PERMUTATION null (2026-09-02). The previous
                    # binomial test treated features as independent Bernoulli
                    # draws, which they emphatically are not -- that is what the
                    # ordination measures -- and it returned p = 0.0 for 10 of
                    # 12 contrasts, the classic tell. Enumerating every distinct
                    # relabelling of the same samples keeps the between-feature
                    # correlation intact, because each permuted split sees the
                    # whole correlated matrix at once.
                    #
                    # Also note the analytic 2/C(n,k) expectation is
                    # CONSERVATIVE on this data: ties from zero-inflation break
                    # separation, so the permutation median runs well below it
                    # and the enrichment ratios are understated, not inflated.
                    # The enumeration INCLUDES the observed labelling, so the
                    # exact conditional p is the plain proportion at least as
                    # extreme -- adding one would double-count the observation
                    # and inflate the floor from 1/126 to 2/127 (the same
                    # error corrected in mwu_exact.mwu_permutation).
                    null = permutation_separations(mat, len(ia))
                    perm_p = float((null >= obs).sum()) / len(null)
                    rows.append({
                        "species": short, "matrix": matrix,
                        "stage_a": sa, "stage_b": sb,
                        "n_a": len(ia), "n_b": len(ib), "n_tested": m,
                        "min_attainable_p": p_min,
                        "k_needed_for_BH_q05": int(np.ceil(p_min * m / 0.05)),
                        "complete_separations": obs,
                        "expected_by_chance": round(exp, 1),
                        "enrichment": round(obs / exp, 2) if exp else np.nan,
                        "perm_null_median": float(np.median(null)),
                        "perm_null_max": int(null.max()),
                        "n_perm_splits": len(null),
                        "perm_p": perm_p,
                        "bh_can_call_any": obs >= int(np.ceil(p_min * m / 0.05)),
                    })
                    print(f"[{short}/{matrix}] {sa} vs {sb}: {obs} separations of {m} "
                          f"(exp {exp:.0f}, {obs/exp:.1f}x, BH needs "
                          f"{int(np.ceil(p_min*m/0.05))}, perm_p={perm_p:.3f})",
                          file=sys.stderr)

    df = pd.DataFrame(rows)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_DIR / "separation_enrichment.tsv", sep="\t", index=False)

    # Figure: enrichment per contrast, with the BH-callability threshold marked.
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, matrix in zip(axes, ["liq", "spore"]):
        sub = df[df["matrix"] == matrix].copy()
        sub["label"] = sub["stage_a"].str[:4] + "-" + sub["stage_b"].str[:4]
        for sp_short, colour in [("dendrobatidis", "#D55E00"),
                                 ("salamandrivorans", "#0072B2")]:
            s = sub[sub["species"] == sp_short]
            x = np.arange(len(s)) + (0.2 if sp_short == "salamandrivorans" else -0.2)
            bars = ax.bar(x, s["enrichment"], width=0.38, color=colour,
                          label=sp_short, edgecolor="white", linewidth=0.5)
            for xi, (_, r) in zip(x, s.iterrows()):
                mark = "*" if r["bh_can_call_any"] else ""
                ax.text(xi, r["enrichment"] + 0.6, f"{r['enrichment']:.0f}x{mark}",
                        ha="center", fontsize=7)
        ax.axhline(1.0, color="#888888", linestyle="--", linewidth=1)
        ax.set_xticks(np.arange(len(s)))
        ax.set_xticklabels(s["label"].tolist(), fontsize=9)
        ax.set_title(f"{matrix} fraction")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("complete-separation enrichment\n(observed / analytic null)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Ordered signal per stage pair, threshold-free (n=5 vs 5).\n"
                 "* = enough simultaneous separations for BH to call any feature at q<0.05",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "separation_enrichment.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT_DIR / "separation_enrichment.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {len(df)} contrasts to {OUT_DIR}/separation_enrichment.tsv", file=sys.stderr)


if __name__ == "__main__":
    main()
