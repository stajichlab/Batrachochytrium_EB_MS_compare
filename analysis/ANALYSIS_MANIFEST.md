# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

## sirius-annotation

- **Status**: transfer + native **full run merged** — 4,957 local features annotated (1,773 transferred + 3,184 native `native-full-e9838293-bagel`; 4,268 with structure; 0 new merge-conflicts from the native fold-in; the pre-existing 69 merged-conflict formula-disagreeing features remain flagged). Transfer: 2,182/2,860 EB SIRIUS-annotated features matched (76.3%). Native pipeline: pilot job `27605104` (149 spectra, 5 shards, completed 2026-08-20, superseded) validated shard size 30; full run job `27718540` (3,927 charge-1+ targets, 3,810 usable spectra, 127 shards) completed 2026-08-25 (one shard, `_11`, TIMEOUT'd on the array's 1h limit and was resubmitted standalone with `--time=02:00:00`, job `27748769_11`), merged via `merged_full/` and folded in 2026-08-25.
- **Purpose**: Attach SIRIUS 6.3.12 formula/structure/CANOPUS annotations to the GNPS2 Everything-Bagel features by m/z+RT+MS2-cosine transfer from the sibling EB project; native SIRIUS run fills the remaining un-annotated charge-1+ features (3,927 targets, ~30 spectra/shard for failure/restart handling) and folds in via `--native-merged`.
- **Key inputs**: `../EB/analysis/sirius_annotation/{sirius_annotations.tsv,sirius_targets.csv,sirius_targets.mgf}`, `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/aligned_features_filled.mgf`, `.../feature_finding_results/aligned_features.csv`.
- **Key outputs**: `analysis/sirius_annotation/{sirius_annotations.tsv,sirius_transfer_map.tsv,sirius_native_targets.csv,sirius_native_targets.mgf,shards_native/,sirius_native_results/,scripts/{import_sirius_transfer.py,select_native_targets.py,export_native_mgf.py,run_sirius_native.sh,run_sirius_native.sbatch}}`.
- **Readme**: [SIRIUS_ANNOTATION.md](sirius_annotation/SIRIUS_ANNOTATION.md)

## ordination

- **Status**: complete (bagel feature table) — all-samples + per-species Bray-Curtis PCoA outputs present (`figures/pcoa_all.{pdf,png}`, `figures/pcoa_condition.{pdf,png}`, `figures/pcoa_stagegroup.{pdf,png}`, `figures/pcoa_axes_all.csv`, `figures/by_species/pcoa_{dendrobatidis,salamandrivorans}.{png,pdf}`); replicates EB F-001 (matrix dominance: all axis1 62.9%, Bd 75.9%, Bsal 70.6% of positive-eigenvalue variance). Pipeline is fully reproducible via `pixi run pcoa-ordination`. New figures color by the sampled condition (6-state `matrix_life_stage`) and by the collapsed life-stage `stage_group` (Zoospore vs Developed).
- **Purpose**: Bray-Curtis PCoA of the feature table across species x matrix (liq/spore) x life_stage (Zoospore/Sporangium/Mature) (GOALS.md goal 3), the Everything-Bagel port of the EB/Rhodotorula pattern.
- **Key inputs**: `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/feature_finding_results/aligned_features.csv`, `data/metdata/curated_gnps_metadata.tsv` → `analysis/ordination/linked_data/` (built by `scripts/build_ordination_table.py`).
- **Key outputs**: `analysis/ordination/{linked_data/{sample_metadata.csv,feature_abundance.csv.gz},figures/{pcoa_axes_all.csv,pcoa_all.png,pcoa_all.pdf,pcoa_condition.png,pcoa_stagegroup.png,by_species/*}}`.
- **Readme**: [ORDINATION.md](ordination/ORDINATION.md)

## differential_features

- **Status**: complete (bagel feature table) — all 30 pairwise contrasts (15 per species) with differential_features.csv.gz + volcano + top_features figures + `comparison_summary.csv`; rank-order concordant with EB (Spearman rho 0.985 Bd / 0.996 Bsal); SIRIUS identity join not yet applied (EB pattern)
- **Purpose**: Feature-level pairwise differential abundance between every (matrix, life_stage) state pair, within each species (GOALS.md goal 2), the Everything-Bagel port of the sibling projects' scripts.
- **Key inputs**: `analysis/ordination/linked_data/` tables (built from the bagel `aligned_features.csv` + curated metadata), `analysis/sirius_annotation/sirius_annotations.tsv` (not yet joined).
- **Key outputs**: `analysis/differential_features/{comparison_summary.csv,<species>_<condA>_vs_<condB>/{differential_features.csv.gz,volcano.png/.pdf,top_features.png/.pdf,top_features.tsv}}`. Repro: `pixi run differential-features`.

## differential_features_primary

- **Status**: complete — 8 primary contrasts with volcano/top_features + SIRIUS-annotated significant-features tables. 4 collapsed life-stage contrasts (Zoospore vs Developed, within species x matrix: liq n=10 vs 20, spore 5 vs 10 per species) + 4 secreted-vs-cellular contrasts (liq vs spore, within species x stage_group). SIRIUS annotation join applied to the merged annotation table (see sirius-annotation, now the full native run); significant feature-rows flagged for liq-vs-spore enrichment (secreted candidates) and curated bioactivity keywords → 103,637 significant rows, 2,639 bioactivity-flagged (re-run 2026-08-25 with the full native SIRIUS annotations).
- **Purpose**: Hypothesis tier on the collapsed life-stage vocabulary (F-002: matrix stratification is required; Sporangium+Mature collapse is power-motivated). Answers "zoospore vs the rest" (life stage) and "secreted fraction vs cell-associated" (matrix) questions behind GOALS.md goals 2-4.
- **Key inputs**: `analysis/ordination/linked_data/` tables, `analysis/sirius_annotation/sirius_annotations.tsv`.
- **Key outputs**: `analysis/differential_features_primary/{primary_comparison_summary.tsv,significant_annotated.tsv,significant_bioactive.tsv,feature_tables_index.html,all_significant_features_summary.tsv,all_significant_features_summary_<species>.html,...}` — note the interactive rollup HTML is chunked per species (dendrobatidis/salamandrivorans) to stay under GitHub's 100 MB per-file limit. Repro: `pixi run differential-features-primary` then `pixi run feature-tables-primary` (sortable/filterable feature-table HTML per contrast + per-species rollups, Rhodotorula strategy).
- **Readme**: [DIFFERENTIAL_FEATURES_PRIMARY.md](differential_features_primary/DIFFERENTIAL_FEATURES_PRIMARY.md)
