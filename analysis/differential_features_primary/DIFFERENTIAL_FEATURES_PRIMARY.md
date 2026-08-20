# DIFFERENTIAL_FEATURES_PRIMARY

## Purpose
Hypothesis tier of the differential analysis, built on the collapsed
life-stage vocabulary established by
`analysis/ordination/scripts/build_ordination_table.py`:

    stage_group = Zoospore | Developed (Sporangium + Mature collapsed)

The 30-way per-condition scan in `analysis/differential_features/` is the
exploratory tier; this is the hypothesis tier, answering the two questions
behind GOALS.md:

1. **Life stage** (within species x matrix): Zoospore vs Developed — "are
   there differences in the features between the zoospore and the rest of
   the stages?". Mature and Sporangium were near-identical in the scan (0
   significant features in every within-matrix stage pair), so collapsing
   them buys power: liq n=10 vs 20, spore n=5 vs 10 per species.
2. **Secreted vs cellular** (within species x stage_group): liq vs spore —
   the liquid-culture supernatant (extracellular/secreted fraction) against
   the spore-pellet (cell-associated) — to isolate features enriched in the
   secreted fraction.

Species are never merged (media is a species confounder: Bd 1% Tryptone vs
Bsal 50% TGHL). All contrasts are stratified by matrix (F-002: matrix
dominates separation and must not be crossed).

## Method (per feature, per comparison)
Identical to the scan tier (`analysis/differential_features/`): prevalence
filter (>=10%) + TSS, Mann-Whitney U, pseudocount-scaled log2 median
fold-change, BH-FDR (direct implementation; no `statsmodels`).

On top of the raw statistics, significant features (q < --fdr) are:
- joined to `analysis/sirius_annotation/sirius_annotations.tsv` (NPC
  pathway/class, ClassyFire class, structure name/SMILES, match status);
- flagged `liq_over_spore_log2fc` = stage-confounded log2(median liq /
  median spore) across that species' whole table — a **secreted-candidate
  hint** (`is_secreted_candidate` when >= +1);
- flagged `bioactive` by a curated keyword regex over structure name + NPC
  pathway/class + ClassyFire class (antibiotic/antimicrobial/mycotoxin/
  siderophore/alkaloid/terpenoid/polyketide/... — a heuristic filter for
  manual curation, not a claim).

## Outputs
`primary_comparison_summary.tsv` (all 8 contrasts, n per side, n_tested,
n_significant) plus one dir per contrast:
`<species>_<groupA>_vs_<groupB>/{differential_features.csv.gz,
volcano.png/.pdf, top_features.png/.pdf, top_features.tsv}`.

Global tables:
- `significant_annotated.tsv` — every significant feature-row across all 8
  contrasts with SIRIUS annotation + secreted-candidate + bioactivity flags.
- `significant_bioactive.tsv` — the bioactivity-flagged subset.

Run: `pixi run differential-features-primary`

## Results (2026-08-20)
8 contrasts, n = 25k-34k features tested each.

| family | contrast (species, groups) | n (a vs b) | n significant |
|---|---|---|---|
| life_stage | Bd liq Zoospore vs Developed | 10 vs 20 | 536 |
| life_stage | Bd spore Zoospore vs Developed | 5 vs 10 | 5,638 |
| life_stage | Bsal liq Zoospore vs Developed | 10 vs 20 | 54 |
| life_stage | Bsal spore Zoospore vs Developed | 5 vs 10 | 7,211 |
| secreted_vs_cellular | Bd liq vs spore (Zoospore) | 10 vs 5 | 17,642 |
| secreted_vs_cellular | Bd liq vs spore (Developed) | 20 vs 10 | 25,051 |
| secreted_vs_cellular | Bsal liq vs spore (Zoospore) | 10 vs 5 | 19,575 |
| secreted_vs_cellular | Bsal liq vs spore (Developed) | 20 vs 10 | 27,930 |

Reads:
- The secreted/cellular matrix is far larger than the life-stage signal in
  the liquid fraction: within liq, only 54-536 features distinguish
  Zoospore from Developed, whereas in the spore fraction 5.6k-7.2k do —
  consistent with F-002 (matrix dominates; the supernatant is compositionally
  dominated by shared media/secreted chemistry).
- 103,637 significant feature-rows cover 33,066 unique features; ~6,688 have
  an NPC class annotation; 3,003 flagged as secreted candidates, 1,041 as
  bioactivity-flagged.
- Bsal (salamandrivorans), which secretes rapidly in culture, shows the
  spore-fraction zoospore-to-developed signal skew (7,211 vs 54 in liq) most
  strongly.

## Caveats
- The bioactivity keyword regex is a filter, not an annotation of record;
  confirm hits against the underlying MS/MS (see `sirius_annotation/`).
- `liq_over_spore_log2fc` is stage-confounded by construction (computed over
  all samples of the species, not within a stage) — use it as a hint, not a
  test; the `secreted_vs_cellular` family is the rigorous test.

## Key inputs
- `analysis/ordination/linked_data/{sample_metadata.csv, feature_abundance.csv.gz}`
- `analysis/sirius_annotation/sirius_annotations.tsv`
