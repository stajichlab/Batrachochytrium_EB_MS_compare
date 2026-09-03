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

### 2026-08-19 — Transfer SIRIUS annotations from the EB project instead of re-running SIRIUS

- **Context**: Ba/Bsal Everything-Bagel features (38,547) need formula/structure/compound-class annotations. A completed SIRIUS 6.3.12 run already exists in the sibling EB project on the MZmine3 feature table sharing the same MassIVE deposit.
- **Decision**: Do not re-run SIRIUS now; transfer the 2,860 EB annotations via m/z + RT + MS2-cosine matching into a project-owned `sirius_annotations.tsv`, and hard-code a `--native-merged` merge path so a future native SIRIUS run upgrades rows in place.
- **Alternatives considered**: (1) Re-run SIRIUS on the Everything-Bagel features now (most correct, but expensive and the user wanted an interim annotation set); (2) try to reuse EB feature IDs directly (impossible — different feature finders, ids do not correspond); (3) m/z+RT-only transfer (insufficient — many isobars within ppm/RT window).
- **Rationale**: The transfer covers 76% of EB's annotated features and is < 15 min to run; native SIRIUS remains the end-state and can overwrite transferred rows since native > transferred by precedence.
- **Consequences**: 1,773 local features annotated interim; 69 `merged_conflict` features require manual review; documentation added at `analysis/sirius_annotation/SIRIUS_ANNOTATION.md`.
**Tags**: metabolomics, sirius, decision, annotation

### 2026-08-19 — Keep one row per local feature id in the accumulated annotation table

- **Context**: The transfer maps multiple EB features onto single local features (296 multi-hit, 69 formula-conflicting). Generating one row per (local feature, EB hit) would fan out joins and double-count abundances.
- **Decision**: `sirius_annotations.tsv` keeps exactly one row per `row ID` (the highest-priority hit by native>transferred, structure-hit, confidence), while `sirius_transfer_map.tsv` audits every EB hit; conflicts are flagged (`merged_conflict`, `n_sirius_hits`, `sirius_hit_ids/formulas`).
- **Alternatives considered**: Wide multi-hit rows; keeping all hits as separate rows.
- **Rationale**: A joinable-by-id table is the cleanest contract for downstream ordination/differential analysis; no information is destroyed because the map retains the full detail.
- **Consequences**: Downstream users must respect `merged_conflict`; the map must not be deleted.
**Tags**: metabolomics, data-model, decision, annotation

### 2026-08-19 — Native SIRIUS run: charge-1+ targets only, small per-shard jobs, pilot-first

- **Context**: The transfer left 4,680 un-annotated features with MS2. User worried the full set was too many; the practical reducers differ from intuition.
- **Decision**: Target the **3,927 singly-charged** (`charge == 1`) un-annotated features (dropping the 753 multi-charged that SIRIUS 6.3.12 cannot process anyway); split into shards of ~30 spectra for fine-grained failure handling/restarts; run a **pilot first** (150 targets, 149 usable, 5 shards) to benchmark per-feature runtime before scaling; keep the array strictly serial (%1) since SIRIUS 6.3.12 serializes login tokens.
- **Alternatives considered**: all 4,680 features (would submit unrunnable multi-charged spectra); 1+ with [M+H]+ only (3,316 — cleaner but 611 fewer); sample-presence filtering (∅, effectively a no-op on this dataset).
- **Rationale**: Charge is the only meaningful reducer here; per-shard independent output project spaces make failures re-submittable without touching other shards; the pilot validates merge/import before the ~3.9k-feature full run.
- **Consequences**: Scripts `select_native_targets.py` / `export_native_mgf.py` / `run_sirius_native.sh`(+`.sbatch`) added; degenerate MGF blocks are dropped with a report. Cost is bounded and restartable.
**Tags**: metabolomics, sirius, decision, native-run, sharding

### 2026-08-20 — Collapse Sporangium+Mature into a "Developed" stage_group for the primary analysis tier

> **CORRECTION (2026-09-02)**: the premise below — "every within-matrix Sporangium-vs-Mature contrast is 0-significant" — is wrong for Bd spore (5,507 significant). The collapse remains justified for Bsal but not for Bd; see the amended F-003. Decision retained as the historical record.

**Context**: The user's main question is Zoospore vs the later stages. The 30-way scan showed every within-matrix Sporangium-vs-Mature pairwise contrast to be 0-significant (F-002), and the user explicitly noted sporangia and mature look very similar — a collapse into two groups buys statistical power (liq 10 vs 20, spore 5 vs 10 per species) for the primary hypothesis tier.

**Decision**: Add `stage_group` (Zoospore | Developed) and `condition_group` (matrix_stage_group) to `linked_data/sample_metadata.csv` in `build_ordination_table.py`; keep the raw 6-state `condition` untouched for the exploratory scan. Add a new `analysis/differential_features_primary/` tier with two contrast families, both stratified by species and matrix: `life_stage` (Zoospore vs Developed, 4 contrasts) and `secreted_vs_cellular` (liq vs spore within stage_group, 4 contrasts).

**Alternatives considered**: Replace the 30-way scan with collapsed contrasts only (rejected — scan stays as exploratory background); collapse via `.str.contains()` hacks in each downstream script (rejected — non-reproducible, scattered).

**Rationale**: Power, reproducibility, and F-002's mandate to stratify by matrix; keeps both the exploratory and hypothesis tiers available.

**Consequences**: Life-stage signal is now read per-fraction: spore fraction carries 5.6k/7.2k FDR-significant features vs 54–536 in liq (F-003); secreted/bioactive targeting uses the liq-vs-spore family plus `is_secreted_candidate` flags on life-stage hits.

**Tags**: analytical-design, life-stage, collapse, differential-abundance, decision, power

### 2026-08-20 — Flag significant features with SIRIUS annotations + curated bioactivity keywords, kept as separate joined tables

**Context**: Downstream goal is to find stage-distinguishing features suggesting secreted or high-bioactivity compounds. The scan tier did not apply the SIRIUS identity join yet (EB pattern).

**Decision**: In the primary tier, join `sirius_annotations.tsv` onto all significant (q<fdr) features across all 8 contrasts into `significant_annotated.tsv`, and write a bioactivity-flagged subset `significant_bioactive.tsv`. Flagging uses a curated keyword regex over structure name + NPC pathway/class + ClassyFire class (antibiotic/antimicrobial/mycotoxin/siderophore/alkaloid/terpenoid/polyketide/...), plus an `is_secreted_candidate` hint = stage-confounded `liq_over_spore_log2fc >= +1`.

**Alternatives considered**: Rely on annotating only volcano-top features (rejected — under-delivers for the goal); mark bioactivity by database lookups only (rejected — no project DB; heuristic regex is the practical interim filter).

**Rationale**: Gives a joinable record of every significant feature with its annotation and context flags without hiding or deleting any rows; caveats document the heuristic as a filter for manual curation, not a claim.

**Consequences**: 103,637 significant feature-rows (33,066 unique features); 6,688 NPC-annotated; 3,003 secreted candidates; 1,041 bioactivity-flagged. Users must treat `bioactive` as reviewing-prompt, and `liq_over_spore_log2fc` as a hint (stage-confounded), not a test.

**Tags**: metabolomics, sirius, bioactivity, secreted-compounds, decision, annotation-join

### 2026-08-20 — Defer the full native SIRIUS run while the Herptile project occupies the `short` queue

**Context**: The 149-spectra pilot (job 27605104, 5 shards × 30, `%1`) validated the full native path (112 formula / 99 structure / 99 NPC / 99 ClassyFire rows merged and folded in). The approved full run (3,921 remaining charge-1+ targets → ~131 shards, ~22–28 h serial on `short`) was ready to launch. But `squeue` shows a sibling project actively running SIRIUS on the same partition: Herptile array `27671478 sirius-herptile-full` (172 remaining serial shards, `%1`) + `27671479 sirius-herptile-full-merge` pending on `afterok:27671478_*`. Herptile's own task logs show SIRIUS login-token (`JOB_WATCHER`) transient failures under load (`Request to Server failed! Try again in 16.0s`).

**Decision**: Hold the full-run submission until `27671478`/`27671479` clear the `short` queue, to avoid two serial SIRIUS arrays competing for queue space and the shared SIRIUS login-token server. Combined with the pilot-first policy, keep `%1` serialization. When launched, generate targets fresh (3,921, excluding the 6 already-annotated features `545/601/1905/2729/8042/8793`) and write to a NEW out-dir `sirius_native_results_full/` (never reuse `sirius_native_results/`, whose `shard_000–004` are a subset of the full targets).

**Alternatives considered**: Launch now and accept contention (rejected — mutual slowdown from token serialization + two arrays on one small queue); launch on a different partition (rejected — `short` is the intended queue and other partitions lack SIRIUS/conda provisioning guarantees for this pipeline).

**Rationale**: The full run has no deadline; starting it today gains nothing but queue/token contention risk and mixes two projects' `JOB_WATCHER` retry storms.

**Consequences**: `SIRIUS_ANNOTATION.md` and `ANALYSIS_MANIFEST.md` marked full-run **deferred (2026-08-20)**; exact launch command recorded. No compute wasted; pilot pipeline stays as the validated template. The 6 duplicate features can later be collapsed in favor of their `native-EB97X` rows.

**Tags**: metabolomics, sirius, decision, native-run, queue-management, deferral

### 2026-08-25 — Confirmed genome-bioactivity-linkage remains blocked on BFD PFAM/antiSMASH, independent of DeepTMHMM completion

**Context**: User reported both the full native SIRIUS run and the DeepTMHMM run as "done" and asked to integrate results. DeepTMHMM output (`TMRs.gff3`, `predicted_topologies.3line`) is genuinely complete and clean for both `dendrobatidis` and `salamandrivorans` (no crash, matches the `set -euo pipefail` fail-loudly fix from 63fa2b3). Ran `pixi run python analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py` live to check whether the full linkage pipeline could now run end-to-end.

**Decision**: Do not attempt to run `build_linkage_tables.py` to completion yet — confirmed it still raises `FileNotFoundError` from `paths.find_bfd_output("pfam_hmmscan", "dendrobatidis", ...)` because BFD has not yet produced `pfam_hmmscan`/`antismash_local` JSON output for either `Batrachochytrium` species locustag (`FCC698BD`, `F61BA062`), even though SignalP/PredGPI/DeepTMHMM (Stage 2 secretion inputs) are all present. Proceeded instead with finishing the SIRIUS native merge (see today's learnings entry) since that portion has no such external blocker.

**Alternatives considered**: Building a DeepTMHMM-only secretion-prediction intermediate table now (via `merge_secretion.py` directly) — deferred as premature scope, since `build_linkage_tables.py` already re-reads all inputs fresh and there is no separate persisted "secretion table" artifact in the current pipeline design to pre-stage.

**Consequences**: `GENOME_BIOACTIVITY_LINKAGE.md`'s existing "Blocked on: the BFD ... run" status still accurately describes reality as of 2026-08-25; no doc change needed there. Re-run `build_linkage_tables.py` once BFD PFAM/antiSMASH output appears under `results/function/pfam_hmmscan/` and `genome_annotation/<out>/antismash_local/` for both species.

**Tags**: genome-bioactivity-linkage, deeptmhmm, sirius, bfd, blocker, pipeline-status

### 2026-08-25 — Superseded: ran genome-bioactivity-linkage end-to-end via local fallback instead of waiting on BFD's own PFAM/antiSMASH/SignalP/PredGPI run

**Context**: The morning's decision (above) held off running `build_linkage_tables.py` because BFD's own shared functional-annotation run had not produced `pfam_hmmscan`/`antismash_local` output for either `Batrachochytrium` locustag. Re-checking disk state found this project already has local fallback results in `analysis/genome_bioactivity_linkage/results/{pfam_hmmscan,antismash_ncbi,predgpi,deeptmhmm,rbh}/<species>/` for both species (from the separate `run_pfam.sh`/`run_signalp.sh`/`run_predgpi.sh` work, commit 413c490) — `find_bfd_output`'s fallback search root (`GBL_ROOT/results/<kind>/<species>`) satisfies the same contract BFD's own run would. Only `signalp/salamandrivorans` was missing (its CPU pass, job `27753328_1`, had `TIMEOUT`'d at 92%).

**Decision**: Do not keep waiting on BFD's own run. Fixed the one missing input directly: rewrote `run_signalp.sh` to use the GPU build (`signalp/6.0h-gpu` on `short_gpu --gres=gpu:1 --time=02:00:00`, was CPU-only on `short`), resubmitted only the missing `salamandrivorans` array index (job `27773831`, completed), then ran `pixi run gbl-build-tables` to completion for both species.

**Alternatives considered**: Keep waiting indefinitely for BFD's shared run (rejected — no ETA, and the local fallback is explicitly designed into `find_bfd_output`/`bfd_antismash_json` for exactly this situation); resubmit SignalP on CPU with more time (rejected once GPU availability was confirmed — GPU removes the bottleneck rather than just buying more wall-clock).

**Rationale**: The fallback path exists specifically so this pipeline is not permanently gated on an external, un-scheduled BFD run; once every Stage 1–2 input the driver needs is present under either search root, there is no reason to keep treating the pipeline as blocked.

**Consequences**: `results/{dendrobatidis,salamandrivorans}_candidate_table.tsv` now exist with real content (0 and 2 candidate rows respectively — see `.living/learnings.md` 2026-08-25 entry on the `is_extracellular` gate for why Bd's count is 0 and not a bug). `GENOME_BIOACTIVITY_LINKAGE.md` and `ANALYSIS_MANIFEST.md` updated to reflect completed-run status; `run_signalp.sh` now defaults to the GPU build for any future rerun (e.g. if BFD's own SignalP output should later supersede the local one, `find_bfd_output`'s search-root order already prefers BFD's path with no code change needed).

**Tags**: genome-bioactivity-linkage, signalp, gpu, bfd, unblocked, pipeline-status

### 2026-08-25 — Relax genome-bioactivity-linkage's `is_extracellular` gate from a hard filter to an informational column

**Context**: The completed end-to-end run (previous decision, same day) confirmed the original strict design — a candidate gene must be `is_extracellular` (SignalP/DeepTMHMM/PredGPI-predicted secreted) to appear in the candidate table at all — produced only 0 (dendrobatidis) and 2 (salamandrivorans) candidate rows on real data, because the domain-hit protein sets (NRPS/PKS/P450/terpene synthase) have essentially zero overlap with the extracellular-scored protein sets. This is biologically expected (these enzymes are cytoplasmic; the metabolite, not the enzyme, is exported), but it meant the pipeline's practical output was almost empty. User asked to relax the gate to surface more candidates for now, with documentation.

**Decision**: Changed `link_compounds_to_genes.build_candidate_table` to accept a `require_extracellular` parameter (default `False`, was implicitly `True`/hardcoded). With the default, every PFAM domain-hit protein of a matched family becomes a candidate regardless of predicted localization; `is_extracellular` remains a column on every output row so it can still be used to filter/prioritize downstream. `require_extracellular=True` restores the original strict behavior for anyone who wants secreted-pathway-only candidates. Reran `pixi run gbl-build-tables`: dendrobatidis 0 → 2,634 rows, salamandrivorans 2 → 8,322 rows.

**Alternatives considered**: Add `is_extracellular` as a tie-breaker/second-order tier factor instead of removing it from the gate entirely (rejected — the design already has an explicit tiered/lexicographic ranking (BGC context, RBH cross-ref) documented as deliberately not a composite score; folding a third signal into tier assignment would blur that scheme, whereas keeping `is_extracellular` as a separate always-present column lets it be filtered on independently without changing the tier semantics). Permanently deleting the strict path (rejected — kept as an opt-in parameter since it's the more conservative, defensible-for-publication reading and may be wanted later).

**Rationale**: The user explicitly asked to see more candidates now; the relaxation is scientifically defensible (most secondary-metabolite biosynthesis genes are legitimately cytoplasmic) and reversible via a single parameter, and the strict-vs-relaxed counts are now documented so nobody mistakes the larger candidate list for higher-confidence evidence than it is.

**Consequences**: `results/{dendrobatidis,salamandrivorans}_candidate_table.tsv` now have substantially more rows; **`tier` alone no longer implies a secreted-pathway gene** — check `is_extracellular` per row if that distinction matters for a given compound (of the current 2,634 + 8,322 rows, only 2 total are `is_extracellular=True`). `GENOME_BIOACTIVITY_LINKAGE.md` (Method step 6, Known caveat #5, Status section) and `ANALYSIS_MANIFEST.md` updated with both the strict and relaxed counts so the before/after is traceable.

**Tags**: genome-bioactivity-linkage, is_extracellular, secretion-prediction, analytical-design, decision, candidate-recovery
