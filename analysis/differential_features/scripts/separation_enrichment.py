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
measure of biological effect. It is why `dendrobatidis_spore_Sporangium_vs_
spore_Mature` reported 5,507 significant on the 38,547-feature table and 0 on
the artifact-filtered 25,157-feature table: the underlying separation counts
barely moved, but removing ~22% duplicate isotope-peak tests pushed the
simultaneous-floor count below the BH threshold. Reporting either number as
"the developmental signal" is an artifact of the denominator.

What this script reports instead
--------------------------------
For each contrast, the number of features showing COMPLETE SEPARATION
(every value in one group strictly exceeds every value in the other), which
is exactly the event that attains p_min, compared against its analytic null
expectation:

    expected = n_tested * 2 / C(n1+n2, n1)
    enrichment = observed / expected
    binomial p = P(X >= observed | n_tested, 2/C(n1+n2,n1))

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
from scipy.stats import binomtest

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
                    bt = binomtest(obs, m, p_min, alternative="greater")
                    rows.append({
                        "species": short, "matrix": matrix,
                        "stage_a": sa, "stage_b": sb,
                        "n_a": len(ia), "n_b": len(ib), "n_tested": m,
                        "min_attainable_p": p_min,
                        "k_needed_for_BH_q05": int(np.ceil(p_min * m / 0.05)),
                        "complete_separations": obs,
                        "expected_by_chance": round(exp, 1),
                        "enrichment": round(obs / exp, 2) if exp else np.nan,
                        "binom_p": bt.pvalue,
                        "bh_can_call_any": obs >= int(np.ceil(p_min * m / 0.05)),
                    })
                    print(f"[{short}/{matrix}] {sa} vs {sb}: {obs} separations of {m} "
                          f"(exp {exp:.0f}, {obs/exp:.1f}x, BH needs "
                          f"{int(np.ceil(p_min*m/0.05))})", file=sys.stderr)

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
