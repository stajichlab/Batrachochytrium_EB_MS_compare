# Codex Guidance

<!-- MYCELIUM:BEGIN -->
This repository uses Mycelium. Before substantial work, read `MYCELIUM.md` and
the routing summary in `.living/INDEX.md`. Follow the post-action protocol in
`MYCELIUM.md` after analysis, data processing, or implementation work.
<!-- MYCELIUM:END -->

# Project Context

Untargeted metabolomics (LC-MS/MS, Q Exactive, positive mode, RP-C18) of the
chytrid fungi *Batrachochytrium dendrobatidis* (Bd) and *Batrachochytrium
salamandrivorans* (Bsal). This repo is a verified analysis workspace built
around a **processed GNPS2 "Everything Bagel" result bundle** (FBMN molecular
networking), NOT around raw spectra.

- **GNPS2 task** `e983829350de4bb39f278cbf22553247`
  (`everything_bagel_workflow`, mode `fbmn`), packaged 2026-08-19.
- **Input/MassIVE deposit** `MSV000090464` (sibling `../MSV000090464`, and
  `../Js_Bd/` for the expanded file-tree export).
- **Bundle root** `data/raw/gnps2_e9838293_bagel/` — read
  `README_FOR_CLAUDE.md` there FIRST. It documents the join key (feature /
  row id), every result table, the annotation workflow and evidence hierarchy,
  artifact screens (PEG, adducts, chimeric spectra), and the metabolomics-USI
  resolver links that render any node's MS/MS live from this task.

Key facts from the curated sample table `data/metdata/curated_gnps_metadata.tsv`:
- 99 biological samples (54 Bd: plates A–F; 45 Bsal: plates V–Z) + 24 QC/IS.
- Per sample: `<plate>_<n>_<type>`, `type` = `liq` | `spore`; a `*C_liq`
  companion (seeding conc 0) is the media control
  (`is_media_control`, `has_C_companion`/`is_C_companion`).
- Three timepoints 8 / 48 / 96 h → `life_stage` Zoospore / Sporangium / Mature.
- `use_in_analysis = True` for 90 rows (the 9 B-plate rows used conditioned
  media — "DO NOT USE"; QC/IS are excluded). Media is a species confounder:
  Bd `1% Tryptone` vs Bsal `50% TGHL` — keep as data, never merge.
- GNPS2 workflow params: `pm_tolerance 0.05`, `fragment_tolerance 0.05`,
  `library_min_matched_peaks 6`, `library_min_similarity 0.7`,
  `detection_preset rare`, `filter_precursor 1`, `topk 1`,
  `formula_prediction_method buddy_fiddle`.

Goals: (1) summarize/characterize the molecular-network components and their
library annotations; (2) tier and curate node annotations (anchor-and-
propagate, mass ladders, class-diagnostic ions); (3) cross the feature/quant
table with the curated metadata for phase/matrix/species comparisons
(timepoint `/zone 8/48/96, liq vs spore, Bd vs Bsal, media controls);
(4) probe molecular families via DeltaMZ cosine-network edges. `data/raw/` is
immutable; derive every processed table reproducibly. See `GOALS.md`.

Reference implementation patterns for differential/volcano dashboards and
ordinations live in the sibling Rhodotorula project
`/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/`
and the Bd-massspec sibling `../EB/`.
