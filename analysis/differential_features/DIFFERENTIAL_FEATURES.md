# DIFFERENTIAL_FEATURES

## Purpose
Pairwise differential-abundance testing of MS2 features between the 6
(matrix, life_stage) states, within each species (15 pairs x 2 species = 30
comparisons) — the Everything-Bagel port of the sibling EB script. Answers
"which features drive the matrix/life-stage separation" behind the ordination
(see `analysis/ordination/ORDINATION.md` and EB F-001).

## Method (per feature, per comparison)
1. Prevalence filter (>=10% of samples across the two groups) + TSS scaling.
2. Mann-Whitney U (rank-based; robust to n=5-10/group, right-skewed areas).
3. Effect size = log2 fold-change of group medians with a feature-scaled
   pseudocount (half the smallest nonzero TSS proportion across both groups).
4. Benjamini-Hochberg FDR across all tested features (implemented directly —
   no `statsmodels` in this env; same order-statistic formula).

## Outputs
`comparison_summary.csv` plus one dir per comparison:
`<species>_<condA>_vs_<condB>/{differential_features.csv.gz, volcano.png/.pdf,
top_features.png/.pdf, top_features.tsv}`.

Run: `pixi run differential-features`

## Results (2026-08-19)
30 comparisons, ~25k-33k features tested each (prevalence >= 10% of the
bigger 38,547-feature table). Matrix contrasts dominate: 16-24k
FDR-significant features per liq-vs-spore pair; liq-vs-liq life-stage
contrasts are ~0 (only Bd liq_Zoospore-vs-liq_Sporangium = 398). Rank-order
concordant with EB (Spearman rho 0.985 Bd / 0.996 Bsal). See
`.living/findings/matrix-dominates-bagel-metabolome.md`.

## Key inputs
- `analysis/ordination/linked_data/{sample_metadata.csv, feature_abundance.csv.gz}`
- `analysis/sirius_annotation/sirius_annotations.tsv` (identity join not yet applied)
