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

### Feature-table HTML (investigation) — `pixi run feature-tables-primary`
Port of the Rhodotorula `generate_compound_table_html.py` strategy
(`analysis/differential_features_primary/scripts/generate_feature_tables.py`):
self-contained, offline-viewable, sortable/filterable tables over the
significant features. Data embedded as JSON; no CDN/no server needed.

- `<comparison>/compound_summary.tsv` + `compound_summary.html` — one per
  contrast, significant rows only (q < 0.05), SIRIUS columns joined.
- `all_significant_features_summary.tsv` — rollup of every significant row
  across all 8 contrasts with `comparison` + `species` columns.
- `all_significant_features_summary_<species>.html` — rollup split into one
  sortable/filterable view per species (`dendrobatidis` 48,867 rows / 59 MB;
  `salamandrivorans` 54,770 rows / 66 MB). Chunked by species so each file
  stays under GitHub's 100 MB per-file hard limit (a single combined rollup
  would be ~125 MB and reject `git push`).
- `feature_tables_index.html` — navigation hub linking all 9 tables.

Design (Rhodotorula-informed): numeric columns sort numerically (not
lexicographic — q-values in sci notation, negative log2FC); identity glyph
(leftmost column, Okabe-Ito): `◇` SIRIUS structure, `○` SIRIUS formula only,
`—` unidentified; ~10 main columns, the ~25 SIRIUS/CANOPUS/transfer-detail
columns in a click-to-expand row panel; filters in priority order = click
sort, q-value `<=` / `|log2FC| >=` numeric thresholds, identity-source chips,
`bioact`/`secreted` checkboxes (the curated flags), free-text search.

Caveat: the per-species rollups embed ~49k–55k rows (~59–66 MB HTML) and the
two largest per-contrast tables ~25–28k rows; they open fine but render slower
than Rhodotorula's ~6k-row tables — prefer filtering or a per-contrast table
when browsing. Regenerate after re-running `differential-features-primary` (the
HTMLs are stale otherwise).

Run: `pixi run feature-tables-primary`

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
- 103,637 significant feature-rows cover 33,066 unique features; 15,791 have
  an NPC class annotation (full native SIRIUS run, 2026-08-25); 62,771
  flagged as secreted candidates, 2,639 as bioactivity-flagged.
- Bsal (salamandrivorans), which secretes rapidly in culture, shows the
  spore-fraction zoospore-to-developed signal skew (7,211 vs 54 in liq) most
  strongly.

## Figures

Volcano (log2FC vs -log10 q) and top-features plots for each of the 8 primary
contrasts (full-resolution `.pdf` alongside each `.png`):

### Life stage (Zoospore vs Developed)

| contrast | volcano | top features |
|---|---|---|
| Bd liq | ![](dendrobatidis_liq_Zoospore_vs_liq_Developed/volcano.png) | ![](dendrobatidis_liq_Zoospore_vs_liq_Developed/top_features.png) |
| Bd spore | ![](dendrobatidis_spore_Zoospore_vs_spore_Developed/volcano.png) | ![](dendrobatidis_spore_Zoospore_vs_spore_Developed/top_features.png) |
| Bsal liq | ![](salamandrivorans_liq_Zoospore_vs_liq_Developed/volcano.png) | ![](salamandrivorans_liq_Zoospore_vs_liq_Developed/top_features.png) |
| Bsal spore | ![](salamandrivorans_spore_Zoospore_vs_spore_Developed/volcano.png) | ![](salamandrivorans_spore_Zoospore_vs_spore_Developed/top_features.png) |

### Secreted vs cellular (liq vs spore)

| contrast | volcano | top features |
|---|---|---|
| Bd Zoospore | ![](dendrobatidis_liq_Zoospore_vs_spore_Zoospore/volcano.png) | ![](dendrobatidis_liq_Zoospore_vs_spore_Zoospore/top_features.png) |
| Bd Developed | ![](dendrobatidis_liq_Developed_vs_spore_Developed/volcano.png) | ![](dendrobatidis_liq_Developed_vs_spore_Developed/top_features.png) |
| Bsal Zoospore | ![](salamandrivorans_liq_Zoospore_vs_spore_Zoospore/volcano.png) | ![](salamandrivorans_liq_Zoospore_vs_spore_Zoospore/top_features.png) |
| Bsal Developed | ![](salamandrivorans_liq_Developed_vs_spore_Developed/volcano.png) | ![](salamandrivorans_liq_Developed_vs_spore_Developed/top_features.png) |

## Caveats
- The bioactivity keyword regex is a filter, not an annotation of record;
  confirm hits against the underlying MS/MS (see `sirius_annotation/`).
- `liq_over_spore_log2fc` is stage-confounded by construction (computed over
  all samples of the species, not within a stage) — use it as a hint, not a
  test; the `secreted_vs_cellular` family is the rigorous test.

## Key inputs
- `analysis/ordination/linked_data/{sample_metadata.csv, feature_abundance.csv.gz}`
- `analysis/sirius_annotation/sirius_annotations.tsv`
