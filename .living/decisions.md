# Decision Log

Append-only log of non-obvious decisions and their rationale.

**Entry template:** copy from `skills/core/templates/decision-log-entry.md` (includes Context, Decision, Alternatives considered, Rationale, Consequences, Tags fields).

### [2026-08-19] Port EB analysis scripts to the Everything-Bagel feature table with a schema-only adapter

**Context**: Goal is a standalone repo reproducing EB's ordination + pairwise differential/volcano analyses, but on the improved Everything-Bagel FBMN feature table (`aligned_features.csv`, 38,547 features) instead of EB's MZMINE3 table (4,107 features). The two tables have different column schemas but the same 90 analysis samples.

**Decision**: Copy `pcoa_ordination.py` and `differential_features.py` verbatim from EB into `analysis/`, and rewrite only `build_ordination_table.py` as a delta from the EB version (rename `row ID/row m/z/row retention time` → `row_id/mz/rt`; map `<stem>.mzML Peak area` → sample_id stem → `_spore`/`_spores`). `linked_data/` schema (row_id, mz, rt, sample cols) is identical to EB, so both downstream scripts run unchanged. Figures are generated locally from this repo's data.

**Alternatives considered**:
- Reuse EB's `gnps_ms2_quant_full.csv` (keep the old table) — rejected: defeats the "improved feature table" purpose; underlying feature data is shifted by Everything-Bagel processing, so results will differ.
- Generalize all three scripts to accept a format-enum flag — rejected: only the build step differs in schema; extra abstraction with no reuse benefit.

**Rationale**: Keeps downstream analyses byte-identical to EB's proven pattern (Bray-Curtis PCoA w/ prevalence≥10% + TSS + fourth-root; Mann-Whitney U + BH-FDR + median log2FC with scaled pseudocount), isolating the difference to the feature-universe shift, which is exactly what we want to compare.

**Consequences**: Results will not be numerically comparable 1:1 (9.8x more features) but rank-order biologically concordant (F-002). Bsal `_spores` stems are preserved verbatim rather than rewritten, so EB-style `spore` plurals must be remembered when joining to other tables.

**Tags**: porting, everything-bagel, ordination, differential-abundance, reproducibility, standalone-repo

## 2026-08-19 — Transfer SIRIUS annotations from the EB project instead of re-running SIRIUS

- **Context**: Ba/Bsal Everything-Bagel features (38,547) need formula/structure/compound-class annotations. A completed SIRIUS 6.3.12 run already exists in the sibling EB project on the MZmine3 feature table sharing the same MassIVE deposit.
- **Decision**: Do not re-run SIRIUS now; transfer the 2,860 EB annotations via m/z + RT + MS2-cosine matching into a project-owned `sirius_annotations.tsv`, and hard-code a `--native-merged` merge path so a future native SIRIUS run upgrades rows in place.
- **Alternatives considered**: (1) Re-run SIRIUS on the Everything-Bagel features now (most correct, but expensive and the user wanted an interim annotation set); (2) try to reuse EB feature IDs directly (impossible — different feature finders, ids do not correspond); (3) m/z+RT-only transfer (insufficient — many isobars within ppm/RT window).
- **Rationale**: The transfer covers 76% of EB's annotated features and is < 15 min to run; native SIRIUS remains the end-state and can overwrite transferred rows since native > transferred by precedence.
- **Consequences**: 1,773 local features annotated interim; 69 `merged_conflict` features require manual review; documentation added at `analysis/sirius_annotation/SIRIUS_ANNOTATION.md`.
- **Tags**: `metabolomics`, `sirius`, `decision`, `annotation`

## 2026-08-19 — Keep one row per local feature id in the accumulated annotation table

- **Context**: The transfer maps multiple EB features onto single local features (296 multi-hit, 69 formula-conflicting). Generating one row per (local feature, EB hit) would fan out joins and double-count abundances.
- **Decision**: `sirius_annotations.tsv` keeps exactly one row per `row ID` (the highest-priority hit by native>transferred, structure-hit, confidence), while `sirius_transfer_map.tsv` audits every EB hit; conflicts are flagged (`merged_conflict`, `n_sirius_hits`, `sirius_hit_ids/formulas`).
- **Alternatives considered**: Wide multi-hit rows; keeping all hits as separate rows.
- **Rationale**: A joinable-by-id table is the cleanest contract for downstream ordination/differential analysis; no information is destroyed because the map retains the full detail.
- **Consequences**: Downstream users must respect `merged_conflict`; the map must not be deleted.
- **Tags**: `metabolomics`, `data-model`, `decision`, `annotation`

## 2026-08-19 — Native SIRIUS run: charge-1+ targets only, small per-shard jobs, pilot-first

- **Context**: The transfer left 4,680 un-annotated features with MS2. User worried the full set was too many; the practical reducers differ from intuition.
- **Decision**: Target the **3,927 singly-charged** (`charge == 1`) un-annotated features (dropping the 753 multi-charged that SIRIUS 6.3.12 cannot process anyway); split into shards of ~30 spectra for fine-grained failure handling/restarts; run a **pilot first** (150 targets, 149 usable, 5 shards) to benchmark per-feature runtime before scaling; keep the array strictly serial (%1) since SIRIUS 6.3.12 serializes login tokens.
- **Alternatives considered**: all 4,680 features (would submit unrunnable multi-charged spectra); 1+ with [M+H]+ only (3,316 — cleaner but 611 fewer); sample-presence filtering (∅, effectively a no-op on this dataset).
- **Rationale**: Charge is the only meaningful reducer here; per-shard independent output project spaces make failures re-submittable without touching other shards; the pilot validates merge/import before the ~3.9k-feature full run.
- **Consequences**: Scripts `select_native_targets.py` / `export_native_mgf.py` / `run_sirius_native.sh`(+`.sbatch`) added; degenerate MGF blocks are dropped with a report. Cost is bounded and restartable.
- **Tags**: `metabolomics`, `sirius`, `decision`, `native-run`, `sharding`
