# Launch the deferred full native SIRIUS run

| Field | Value |
|-------|-------|
| **Date** | 2026-08-20 |
| **Author** | jstajich (via agent session) |
| **Priority** | medium |
| **Status** | done (2026-08-25) — full run completed, job 27718540 (127 shards, one resubmitted standalone as 27748769_11 after TIMEOUT), merged and folded into `sirius_annotations.tsv`; see `analysis/sirius_annotation/SIRIUS_ANNOTATION.md` |
| **Category** | infrastructure |
| **Related analyses** | `analysis/sirius_annotation/` |
| **Related data** | `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/feature_finding_results/aligned_features.{csv,mgf:filled}` |

## Description

Run native SIRIUS 6.3.12 on the remaining un-annotated charge-1+ features of the Everything-Bagel feature table, reusing the validated pilot/shared pipeline. Deferred (2026-08-20) because the sibling Herptile array `27671478` (+ merge `27671479`) occupies the `short` queue and contends for the shared SIRIUS login-token server.

## Motivation

Annotates the ~3,921 un-annotated singly-charged features with formula/structure/NPC/ClassyFire, upgrading rows in place with `annotation_origin = 'native'` > `transferred` precedence.

## Proposed Approach

1. Wait for `27671478` + `27671479` to exit (`squeue -u $USER -p short`; currently task 7/176 running).
2. Regenerate targets with `select_native_targets.py` (should now exclude the 6 already-annotated features `545/601/1905/2729/8042/8793` → ~3,921 targets).
3. Launch `bash analysis/sirius_annotation/scripts/run_sirius_native.sh 30 0 1` writing to a FRESH out-dir `sirius_native_results_full/` (never reuse `sirius_native_results/` — its `shard_000–004` are a pilot subset of the full targets).
4. Merge shards via EB `merge_sirius_shards.py`, import via `import_sirius_transfer.py --native-merged ... --native-label native-EB97X` (only if not already done by a concurrent agent).
5. Collapse the 6 duplicate-row features in favor of their `native-EB97X` rows.
6. Re-run downstream differential fold-ins; update `SIRIUS_ANNOTATION.md` + `ANALYSIS_MANIFEST.md` with measured runtime.

## Acceptance Criteria

- [ ] Array completes on `short`, all shards COMPLETED (≈131 shards, `%1`, ~22–28 h)
- [ ] Merged + imported; `sirius_annotations.tsv` grows by the native fold-in with 0 new merge-conflicts
- [ ] Doc/manifest updated with actual runtime and row counts

## Notes

- Pilot benchmark: 30-spectra shard = 7:38–13:46 wall, MaxRSS 6.2–12.8 GB → ~15–25 s/feature.
- Keep changes to the shared `sirius_annotations.tsv` minimal/coordinated with concurrent agents; repo is being reshaped by a concurrent git effort.
- Decisions: `.living/decisions.md` 2026-08-20 (deferral); learnings: `.living/learnings.md` 2026-08-20 (queue/login-token contention).
