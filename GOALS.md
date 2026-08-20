# GOALS.md — Batrachochytrium Everything-Bagel Metabolomics

Untargeted LC-MS/MS metabolomics of the chytrid fungi *Batrachochytrium
dendrobatidis* (Bd) and *Batrachochytrium salamandrivorans* (Bsal). This repo
is the analysis workspace for a **processed GNPS2 "Everything Bagel" (FBMN
molecular networking) result bundle** from task
`e983829350de4bb39f278cbf22553247` on MassIVE `MSV000090464`.

## Data at a Glance

- **GNPS2 Everything Bagel bundle** in `data/raw/gnps2_e9838293_bagel/`
  (read `README_FOR_CLAUDE.md` there first):
  - `aligned_features.csv` — feature × sample quantification (the primary
    quant table; feature/row id = join key).
  - `aligned_features_filled.mgf` — one MS/MS per feature (`SCANS` =
    feature id).
  - `merged_feature_library_search_results.tsv` — spectral-library hits,
    keyed on feature id.
  - `filtered_pairs.tsv` — molecular-network edges (`CLUSTERID1/2`,
    `DeltaMZ`, `Cosine`, `MatchedPeaks`).
  - `network.graphml` / `network_singletons.graphml` — consolidated node +
    edge + annotation view (`component` = molecular family).
- **Workflow params**: `pm_tolerance 0.05`, `fragment_tolerance 0.05`,
  `library_min_matched_peaks 6`, `library_min_similarity 0.7`,
  `detection_preset rare`, `filter_precursor 1`, `topk 1`,
  `formula_prediction_method buddy_fiddle`; mode `fbmn`.
- **Curated sample table** `data/metdata/curated_gnps_metadata.tsv`: 123 rows
  = 99 biological samples (54 Bd, 45 Bsal) + 24 QC/IS. Design per sample:
  `<plate>_<n>_<type>` (`liq` | `spore`), a `*C_liq` media-control companion
  (seeding conc 0), timepoints 8/48/96 h →
  Zoospore/Sporangium/Mature, media Bd `1% Tryptone` vs Bsal `50% TGHL`
  (species confounder — keep as data, never merge).

## Goals

1. **Characterize the molecular network.** Summarize the largest components
   (molecular families), their sizes, edge structure (DeltaMZ / cosine ladders),
   and library annotations. List confidently library-annotated nodes
   (high `MQScore`, many `SharedPeaks`) and the compound classes present.
2. **Tier & curate annotations.** Categorize node annotations by evidence
   (library-verified > network-propagated-from-anchor > diagnostic-ion only).
   Run anchor-and-propagate within components, mass-ladder unnamed members,
   and screen the whole dataset for artifacts (PEG/polymer ladders, adduct
   pairs, chimeric/co-isolated spectra, mass-defect outliers) before trusting
   any family-level call.
3. **Sample-level comparisons.** Cross the feature/quant table with the
   curated metadata for: life-stage / timepoint (8/48/96 h − Zoospore /
   Sporangium / Mature), matrix (liq vs spore), species (Bd vs Bsal), and
   media-control subtraction. Use only `use_in_analysis = True` rows;
   TIC-normalize, prefer presence/absence on sparse features.
4. **Probe molecular families.** Use DeltaMZ cosine-network edges to propose
   modifications / analogs within a component (CH2 ladder, H2
   saturation/desaturation, C2H2, conjugations); validate each proposal
   against the fragment evidence and RT trend.
5. **Reproducibility.** All derived tables/figures from committed scripts run
   against the immutable bundle (`data/raw/` never modified). Link spectra via
   the task's metabolomics-USI resolver so every claim is checkable live.

## Non-Goals

- Re-processing raw `.raw`/`.mzML` files (the bundle is the processed truth).
- Re-running GNPS2 feature finding / networking (use this run's outputs
  unless a new submission is explicitly requested).

## Out of Scope / Deferred (see `todo/TODO_REGISTRY.md`)

- Deeper structure annotation (SIRIUS 6.x, ModiFinder) — sibling project
  `../EB/` covers the SIRIUS reprocessing track.
- MS1 isotope-envelope-backed halogen / charge conclusions (not in this
  export — no MS1 envelopes).
