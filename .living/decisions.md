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

## 2026-08-20 — Collapse Sporangium+Mature into a "Developed" stage_group for the primary analysis tier

**Context**: The user's main question is Zoospore vs the later stages. The 30-way scan showed every within-matrix Sporangium-vs-Mature pairwise contrast to be 0-significant (F-002), and the user explicitly noted sporangia and mature look very similar — a collapse into two groups buys statistical power (liq 10 vs 20, spore 5 vs 10 per species) for the primary hypothesis tier.

**Decision**: Add `stage_group` (Zoospore | Developed) and `condition_group` (matrix_stage_group) to `linked_data/sample_metadata.csv` in `build_ordination_table.py`; keep the raw 6-state `condition` untouched for the exploratory scan. Add a new `analysis/differential_features_primary/` tier with two contrast families, both stratified by species and matrix: `life_stage` (Zoospore vs Developed, 4 contrasts) and `secreted_vs_cellular` (liq vs spore within stage_group, 4 contrasts).

**Alternatives considered**: Replace the 30-way scan with collapsed contrasts only (rejected — scan stays as exploratory background); collapse via `.str.contains()` hacks in each downstream script (rejected — non-reproducible, scattered).

**Rationale**: Power, reproducibility, and F-002's mandate to stratify by matrix; keeps both the exploratory and hypothesis tiers available.

**Consequences**: Life-stage signal is now read per-fraction: spore fraction carries 5.6k/7.2k FDR-significant features vs 54–536 in liq (F-003); secreted/bioactive targeting uses the liq-vs-spore family plus `is_secreted_candidate` flags on life-stage hits.

**Tags**: `analytical-design`, `life-stage`, `collapse`, `differential-abundance`, `decision`, `power`

## 2026-08-20 — Flag significant features with SIRIUS annotations + curated bioactivity keywords, kept as separate joined tables

**Context**: Downstream goal is to find stage-distinguishing features suggesting secreted or high-bioactivity compounds. The scan tier did not apply the SIRIUS identity join yet (EB pattern).

**Decision**: In the primary tier, join `sirius_annotations.tsv` onto all significant (q<fdr) features across all 8 contrasts into `significant_annotated.tsv`, and write a bioactivity-flagged subset `significant_bioactive.tsv`. Flagging uses a curated keyword regex over structure name + NPC pathway/class + ClassyFire class (antibiotic/antimicrobial/mycotoxin/siderophore/alkaloid/terpenoid/polyketide/...), plus an `is_secreted_candidate` hint = stage-confounded `liq_over_spore_log2fc >= +1`.

**Alternatives considered**: Rely on annotating only volcano-top features (rejected — under-delivers for the goal); mark bioactivity by database lookups only (rejected — no project DB; heuristic regex is the practical interim filter).

**Rationale**: Gives a joinable record of every significant feature with its annotation and context flags without hiding or deleting any rows; caveats document the heuristic as a filter for manual curation, not a claim.

**Consequences**: 103,637 significant feature-rows (33,066 unique features); 6,688 NPC-annotated; 3,003 secreted candidates; 1,041 bioactivity-flagged. Users must treat `bioactive` as reviewing-prompt, and `liq_over_spore_log2fc` as a hint (stage-confounded), not a test.

**Tags**: `metabolomics`, `sirius`, `bioactivity`, `secreted-compounds`, `decision`, `annotation-join`
