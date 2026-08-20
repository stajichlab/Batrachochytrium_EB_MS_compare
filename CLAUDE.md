# Claude Code Guidance

<!-- MYCELIUM:BEGIN -->
This repository uses Mycelium. Before substantial work, read `MYCELIUM.md` and
the routing summary in `.living/INDEX.md`. Follow the post-action protocol in
`MYCELIUM.md` after analysis, data processing, or implementation work.
<!-- MYCELIUM:END -->

# Project Context

Untargeted metabolomics (LC-MS/MS) of the chytrid fungi *Batrachochytrium
dendrobatidis* (Bd) and *Batrachochytrium salamandrivorans* (Bsal). This repo
analyzes a **processed GNPS2 "Everything Bagel" (FBMN molecular networking)
result bundle** — the downstream truth tables, not raw spectra.

- **GNPS2 task** `e983829350de4bb39f278cbf22553247`
  (`everything_bagel_workflow`, mode `fbmn`), Input **MassIVE `MSV000090464`**.
- **Bundle root** `data/raw/gnps2_e9838293_bagel/` — read
  `README_FOR_CLAUDE.md` there FIRST (join key = feature/row id; every result
  table; annotation + artifact-screen guidance; metabolomics-USI resolver
  links to render any node's MS/MS).
- **Curated sample table** `data/metdata/curated_gnps_metadata.tsv`: 99
  biological (54 Bd plates A–F, 45 Bsal plates V–Z) + 24 QC/IS; timepoints
  8/48/96 h → Zoospore/Sporangium/Mature; `liq`/`spore` matrix with `*C_liq`
  media controls; `use_in_analysis = True` for 90 rows (B plates = conditioned
  media, "DO NOT USE"; QC/IS excluded). Media is a species confounder
  (Bd `1% Tryptone` vs Bsal `50% TGHL`) — keep as data.
- Goals: summarize molecular-network components + library annotations, tier
  node annotations, cross feature quant with the curated metadata for
  phase/matrix/species comparisons, probe families via DeltaMZ edges. See
  `AGENTS.md`, `GOALS.md`, `MYCELIUM.md`.
- Reference implementation patterns (differential/volcano dashboards,
  ordinations, SIRIUS pipelines) live in the sibling Rhodotorula project
  `/bigdata/stajichlab/shared/projects/Rhodotorula/Rhodotorula_Metabolites/Rhodotorula_pheno_MS/`
  and the Bd-massspec sibling `../EB/`.
