#!/usr/bin/env python3
"""Media-consumption analysis: which medium components does each species deplete?

Why
---
The corrected re-analysis (CORRECTED_REANALYSIS_REPORT.md) turned the medium
from a nuisance into the signal. Two independent lines -- residue composition
of the blank-clearing peptides (Pro+Hyp 24.3% Bd / 16.9% Bsal, casein/gelatin
levels) and the molecular-family tier (Bsal's blank-enriched components are
proline-rich peptide families) -- point at secreted proteolysis of medium
protein rather than secondary-metabolite biosynthesis.

That makes the *complement* of the usual analysis directly informative:
features **depleted** relative to their own media blank are medium components
the fungus consumed or transformed. This is arguably the most robust
liquid-fraction readout available here, because it does not require the
feature to be fungal in origin -- only that the fungus changed its abundance.

It also yields a falsifiable prediction from the repo's own comparative
genomics. Bsal carries a large MEROPS **M36 fungalysin** expansion (328 M36
hits vs 39 in Bd; 233/247 secreted protease candidates -- PROTEASE_CANDIDATES.md,
literature-corroborated by Yu et al. 2025). If that expansion is functional in
these cultures, Bsal should deplete peptide-class medium features **more** than
Bd does.

Caveat that bounds the cross-species comparison: the media DIFFER (Bd 1%
tryptone, Bsal 50% TGHL) and were run two months apart, so absolute depletion
counts are not comparable across species. What IS comparable is the
*compositional* question -- within each species, is the depleted set enriched
for peptides relative to that species' own tested background? That is a
within-species contrast and the medium confound cancels.

Method
------
Per species x life_stage, paired by (plate, replicate) exactly as
`background_subtraction` does:

  depleted  : fungal <= (1/min_fc) x its own plate blank in >= 4 of 5 plates
  enriched  : fungal >= min_fc x its own plate blank in >= 4 of 5 plates  (the
              existing blank-clearing rule, recomputed here for symmetry)

Then, within each species, Fisher-exact test whether the depleted set is
enriched for SIRIUS "Amino acids and Peptides" relative to all tested features
of that species. Also reports the stage trend (does depletion grow
Zoospore -> Mature, i.e. progressive consumption?).

Usage: python3 scripts/media_consumption.py
Outputs: analysis/differential_features_primary/media_consumption/
         {consumption.tsv, class_enrichment.tsv, consumption.png}
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

REPO = Path(__file__).resolve().parents[3]
LINKED = REPO / "analysis" / "ordination" / "linked_data"
SIRIUS = REPO / "analysis" / "sirius_annotation" / "sirius_annotations.tsv"
OUT = REPO / "analysis" / "differential_features_primary" / "media_consumption"

sys.path.insert(0, str(REPO / "analysis" / "genome_bioactivity_linkage" / "scripts"))
import background_subtraction as bs  # noqa: E402

SPECIES = ["Batrachochytrium dendrobatidis", "Batrachochytrium salamandrivorans"]
STAGES = ["Zoospore", "Sporangium", "Mature"]
MIN_FC = 2.0
MIN_PLATES = bs._MIN_PAIRED_PLATES
PEPTIDE_CLASS = "Amino acids and Peptides"


def paired_counts(features: pd.DataFrame, pairs: list[tuple[str, str]]):
    """(n_plates_depleted, n_plates_enriched) per feature, using the paired rule."""
    pc = bs._PSEUDOCOUNT
    dep = np.zeros(len(features), dtype=int)
    enr = np.zeros(len(features), dtype=int)
    for f_col, b_col in pairs:
        ratio = (features[f_col] + pc) / (features[b_col] + pc)
        dep += (ratio <= 1.0 / MIN_FC).to_numpy(dtype=int)
        enr += (ratio >= MIN_FC).to_numpy(dtype=int)
    return dep, enr


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    features = bs.load_feature_intensities()
    meta = bs.load_metadata()

    # Restrict to the artifact-filtered analysis-matrix universe.
    kept = set(pd.read_csv(LINKED / "feature_abundance.csv.gz", usecols=["row_id"])["row_id"].astype(int))
    features = features.loc[features.index.isin(kept)]
    print(f"features in analysis-matrix universe: {len(features)}", file=sys.stderr)

    sirius = pd.read_csv(SIRIUS, sep="\t")
    sirius["row ID"] = sirius["row ID"].astype(int)
    npc = dict(zip(sirius["row ID"], sirius["sirius_npc_pathway"]))

    rows, per_species_sets = [], {}
    for sp in SPECIES:
        short = sp.split()[-1]
        dep_union: set[int] = set()
        enr_union: set[int] = set()
        for stage in STAGES:
            scoped = meta[
                (meta["species"] == sp) & (meta["life_stage"] == stage)
                & (meta["matrix"] == "liq") & (meta["use_in_analysis"] == True)  # noqa: E712
            ]
            pairs = bs._plate_pairs(scoped, features)
            if not pairs:
                sys.exit(f"no plate pairs resolved for {short}/{stage}")
            dep, enr = paired_counts(features, pairs)
            thresh = min(MIN_PLATES, len(pairs))
            dep_ids = set(features.index[dep >= thresh].astype(int))
            enr_ids = set(features.index[enr >= thresh].astype(int))
            dep_union |= dep_ids
            enr_union |= enr_ids
            rows.append({
                "species": short, "life_stage": stage, "n_plate_pairs": len(pairs),
                "n_tested": len(features),
                "n_depleted": len(dep_ids), "n_enriched": len(enr_ids),
                "frac_depleted": round(len(dep_ids) / len(features), 4),
                "frac_enriched": round(len(enr_ids) / len(features), 4),
            })
            print(f"[{short}/{stage}] pairs={len(pairs)} depleted={len(dep_ids)} "
                  f"enriched={len(enr_ids)} of {len(features)}", file=sys.stderr)
        per_species_sets[short] = (dep_union, enr_union)

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "consumption.tsv", sep="\t", index=False)

    # Within-species compositional test: is the depleted set peptide-enriched?
    enr_rows = []
    for short, (dep_union, enr_union) in per_species_sets.items():
        annotated = {f for f in features.index.astype(int) if pd.notna(npc.get(f))}
        bg_pep = sum(1 for f in annotated if npc[f] == PEPTIDE_CLASS)
        bg_tot = len(annotated)
        for label, ids in [("depleted", dep_union), ("enriched", enr_union)]:
            ann = ids & annotated
            k = sum(1 for f in ann if npc[f] == PEPTIDE_CLASS)
            n = len(ann)
            if n == 0:
                continue
            odds, p = fisher_exact(
                [[k, n - k], [bg_pep - k, (bg_tot - n) - (bg_pep - k)]],
                alternative="greater")
            enr_rows.append({
                "species": short, "set": label, "n_annotated_in_set": n,
                "n_peptide_in_set": k, "frac_peptide_in_set": round(k / n, 4),
                "background_frac_peptide": round(bg_pep / bg_tot, 4),
                "odds_ratio": round(odds, 3), "fisher_p": p,
            })
            print(f"  {short} {label}: peptide {k}/{n} = {k/n:.1%} vs background "
                  f"{bg_pep/bg_tot:.1%}  OR={odds:.2f} p={p:.2e}", file=sys.stderr)
    ce = pd.DataFrame(enr_rows)
    ce.to_csv(OUT / "class_enrichment.tsv", sep="\t", index=False)

    # Figure: depleted vs enriched counts by stage, per species.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, short in zip(axes, ["dendrobatidis", "salamandrivorans"]):
        s = df[df["species"] == short]
        x = np.arange(len(s))
        ax.bar(x - 0.2, s["n_depleted"], width=0.38, color="#0072B2",
               label="depleted vs blank", edgecolor="white")
        ax.bar(x + 0.2, s["n_enriched"], width=0.38, color="#D55E00",
               label="enriched vs blank", edgecolor="white")
        for xi, (_, r) in zip(x, s.iterrows()):
            ax.text(xi - 0.2, r["n_depleted"] + 30, f"{int(r['n_depleted'])}",
                    ha="center", fontsize=8)
            ax.text(xi + 0.2, r["n_enriched"] + 30, f"{int(r['n_enriched'])}",
                    ha="center", fontsize=8)
        ax.set_xticks(x)
        ax.set_xticklabels(s["life_stage"], fontsize=9)
        ax.set_title(short)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(f"features ({MIN_FC:g}x, >= {MIN_PLATES}/5 plates)")
    axes[0].legend(frameon=False, fontsize=9)
    fig.suptitle("Medium components depleted vs released, per life stage\n"
                 "(paired to each well's own C_liq blank; media differ between "
                 "species so compare stages within a species, not across)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "consumption.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "consumption.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {len(df)} stage rows and {len(ce)} enrichment rows to {OUT}",
          file=sys.stderr)


if __name__ == "__main__":
    main()
