---
session_id: 2026-09-02-001
project: batrachochytrium-eb-ms-compare
branch: "main"
started: 2026-09-02T22:35:19-0700
ended:
duration_minutes:
files_changed:
---

## Session Log

### 22:35 — Session started
- Branch: `main`
- Resuming from: (first session)

### 23:15 — Corrected re-analysis committed; global knowledge initialized
- **Command**: `pixi run build-ordination-table / pcoa-ordination / differential-features / separation-enrichment / differential-features-primary / feature-tables-primary / lifestage-trend / usi-curation / peptide-origin`; `init_knowledge.py`
- **Result**: five defects fixed and full pipeline re-run on 60 fungal samples x 25,157 features (was 90 x 38,547). All 12 stage pairs shown 3.0-42.6x separation-enriched with none BH-callable; blank-clearing MS2-backed liq peptides shown proline-rich at casein/gelatin levels (Bd 24.3%, Bsal 16.9%; proteome baseline rejected p=5.8e-12 / 5.3e-10). F-002 amended, F-003 partially retracted. Committed 8d8f9e6 + 2555ab7. `~/.mycelium/knowledge/` bootstrapped (15 domain files) and two transferable statistics learnings promoted to data-analysis.md.
- **Output**: `CORRECTED_REANALYSIS_REPORT.md`; `analysis/differential_features/{separation_enrichment.*,scripts/mwu_exact.py}`; `analysis/differential_features_primary/{lifestage_trend/,peptide_origin/}`
