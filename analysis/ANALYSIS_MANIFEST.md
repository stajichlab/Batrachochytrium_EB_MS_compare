# Analysis Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

## sirius-annotation

- **Status**: transfer complete; **native run in progress (pilot)** — transfer: 2,182/2,860 EB SIRIUS-annotated features matched (76.3%); 1,773 local features annotated (1,773 formula, 1,657 structure); 69 merged-conflict (formula-disagreeing) features flagged. Native pipeline built (target selection → MGF export → sharding → SLURM array), pilot job `27605104` (149 usable spectra, 5 shards) running to benchmark runtime.
- **Purpose**: Attach SIRIUS 6.3.12 formula/structure/CANOPUS annotations to the GNPS2 Everything-Bagel features by m/z+RT+MS2-cosine transfer from the sibling EB project; native SIRIUS run fills the remaining un-annotated charge-1+ features (3,927 targets, ~30 spectra/shard for failure/restart handling) and folds in via `--native-merged`.
- **Key inputs**: `../EB/analysis/sirius_annotation/{sirius_annotations.tsv,sirius_targets.csv,sirius_targets.mgf}`, `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/aligned_features_filled.mgf`, `.../feature_finding_results/aligned_features.csv`.
- **Key outputs**: `analysis/sirius_annotation/{sirius_annotations.tsv,sirius_transfer_map.tsv,sirius_native_targets.csv,sirius_native_targets.mgf,shards_native/,sirius_native_results/,scripts/{import_sirius_transfer.py,select_native_targets.py,export_native_mgf.py,run_sirius_native.sh,run_sirius_native.sbatch}}`.
- **Readme**: [SIRIUS_ANNOTATION.md](sirius_annotation/SIRIUS_ANNOTATION.md)

## ordination

- **Status**: complete (bagel feature table) — all-samples + per-species Bray-Curtis PCoA outputs present (`figures/pcoa_all.{pdf,png}`, `figures/pcoa_axes_all.csv`, `figures/by_species/pcoa_{dendrobatidis,salamandrivorans}.{png,pdf}`); replicates EB F-001 (matrix dominance: all axis1 62.9%, Bd 75.9%, Bsal 70.6% of positive-eigenvalue variance). Pipeline is fully reproducible via `pixi run pcoa-ordination`.
- **Purpose**: Bray-Curtis PCoA of the feature table across species x matrix (liq/spore) x life_stage (Zoospore/Sporangium/Mature) (GOALS.md goal 3), the Everything-Bagel port of the EB/Rhodotorula pattern.
- **Key inputs**: `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/feature_finding_results/aligned_features.csv`, `data/metdata/curated_gnps_metadata.tsv` → `analysis/ordination/linked_data/` (built by `scripts/build_ordination_table.py`).
- **Key outputs**: `analysis/ordination/{linked_data/{sample_metadata.csv,feature_abundance.csv.gz},figures/{pcoa_axes_all.csv,pcoa_all.png,pcoa_all.pdf,by_species/*}}`.
- **Readme**: (none yet — see `.living/findings/matrix-dominates-bagel-metabolome.md`)

## differential_features

- **Status**: complete (bagel feature table) — all 30 pairwise contrasts (15 per species) with differential_features.csv.gz + volcano + top_features figures + `comparison_summary.csv`; rank-order concordant with EB (Spearman rho 0.985 Bd / 0.996 Bsal); SIRIUS identity join not yet applied (EB pattern)
- **Purpose**: Feature-level pairwise differential abundance between every (matrix, life_stage) state pair, within each species (GOALS.md goal 2), the Everything-Bagel port of the sibling projects' scripts.
- **Key inputs**: `analysis/ordination/linked_data/` tables (built from the bagel `aligned_features.csv` + curated metadata), `analysis/sirius_annotation/sirius_annotations.tsv` (not yet joined).
- **Key outputs**: `analysis/differential_features/{comparison_summary.csv,<species>_<condA>_vs_<condB>/{differential_features.csv.gz,volcano.png/.pdf,top_features.png/.pdf,top_features.tsv}}`. Repro: `pixi run differential-features`.
