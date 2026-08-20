# Data Manifest

<!-- Add entries below using the appropriate manifest entry template. -->

### gnps2-everything-bagel (FBMN molecular network result bundle)
```yaml
name: gnps2-everything-bagel
type: metabolomics
source: GNPS2 task e983829350de4bb39f278cbf22553247 (everything_bagel_workflow, mode fbmn); packaged 2026-08-19; input MassIVE deposit MSV000090464
date_acquired: 2026-08-19
format: CSV (feature finding / library search / quant), MGF (consensus spectra), TSV (network edges), GraphML (networks), YAML (workflow params)
size: 234M (data/raw/gnps2_e9838293_bagel/)
raw_path: data/raw/gnps2_e9838293_bagel/
status: raw
known_issues:
  - Bsal filenames use plural `_spores` in area columns; Bd plates A-F use singular `_spore` (directional species convention, not a typo -- match stems, never rewrite)
  - Feature table contains full 38,547 features incl. low-prevalence / sparse rows; analysis scripts apply prevalence filters downstream
  - data/raw/ is immutable; all derived tables must be regenerated reproducibly from here
access_restrictions: none
tags: [gnps2, fbmn, molecular-network, ms2, massspec, chytrid, bd, bsal]
```

Key bundle paths:
- `README_FOR_CLAUDE.md` — bundle documentation; read first (join key = feature/row id, annotation workflow, artifact screens)
- `nf_output/feature_finding/feature_finding_results/aligned_features.csv` — aligned feature/quant table: 38,547 features x 152 cols (`row ID`, `row m/z`, `row retention time`, ... + 123 `<sample>.mzML Peak area` columns incl. 24 QC/IS)
- `nf_output/feature_finding/feature_finding_results/aligned_features_filled.mgf` — consensus MS/MS spectra (SCANS = feature id)
- `nf_output/feature_library_search/` — merged library search results + GNPS annotations
- `nf_output/networking/` — `filtered_pairs.tsv` / `merged_pairs.tsv` (edges: CLUSTERID1/2, DeltaMZ, Cosine, MatchedPeaks), `network.graphml` (~16M), `network_singletons.graphml` (~74M)

### curated-gnps-metadata
```yaml
name: curated-gnps-metadata
type: metabolomics
source: Manually curated from MassIVE MSV000090464 submission metadata (MSV000090464 Batrachochytrium/curated_gnps_metadata.tsv at submit time)
date_acquired: 2026-08-19
format: TSV
rows: 123
columns: 49
size: small
raw_path: data/raw/gnps2_e9838293_bagel/ (submitted copy)
metadata_path: data/metdata/curated_gnps_metadata.tsv
status: validated
known_issues:
  - 9 B-plate rows used conditioned media (DO NOT USE / use_in_analysis = False); QC + IS rows also excluded from analysis
  - Media is a species confounder: Bd 1% Tryptone vs Bsal 50% TGHL -- keep as data, never merge
  - Directory spelled `metdata` (project-local convention, distinct from EB's `metadata/`)
access_restrictions: none
tags: [metadata, samples, bd, bsal, timepoint, matrix]
```

Notes:
- 99 biological samples (54 Bd plates A-F, 45 Bsal plates V-Z) + 24 QC/IS; timepoints 8/48/96 h (Zoospore/Sporangium/Mature); matrix liq|spore with `*C_liq` media-control companions.
- `use_in_analysis = True` for 90 rows — these match the feature table's non-QC/IS area columns exactly (123/123 stems align).

### bagel-ordination-linked-data (derived)
```yaml
name: bagel-ordination-linked-data
type: metabolomics
source: analysis/ordination/scripts/build_ordination_table.py -- joins curated-gnps-metadata (use_in_analysis==True) to gnps2-everything-bagel aligned_features.csv area columns
date_acquired: 2026-08-19
format: CSV / CSV.GZ
rows: 90 samples (sample_metadata.csv) x 38,547 features (feature_abundance.csv.gz)
status: derived
known_issues:
  - Bsal `_spores` stems preserved verbatim from both sources (consistent, no rewriting)
  - Raw peak area; zero-filled where feature absent
processed_path: analysis/ordination/linked_data/
tags: [derived, linked-data, ordination, differential]
```
