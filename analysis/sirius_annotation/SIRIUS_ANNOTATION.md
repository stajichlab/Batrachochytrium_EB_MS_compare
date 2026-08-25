# SIRIUS Annotation (Transferred + Native Merge Path)

## Purpose

Attach formula / structure / compound-class (CANOPUS, NPC, ClassyFire)
annotations from SIRIUS 6.3.12 to this project's GNPS2 "Everything Bagel"
features (38,547 features, MassIVE `MSV000090464`).

SIRIUS was run on the **sibling EB project**
(`/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/`) on the
MZmine3 feature table (4,107 features). Because the two projects share the
deposit but used **different feature finders**, feature IDs do not correspond
— so annotations are **transferred by matching** precursor `m/z` + RT + MS2
spectral cosine, rather than re-running SIRIUS. This document records how the
transfer works, what it produced, its caveats, and how to fold in the next
**native** SIRIUS run when one is completed on these features.

## Status

**status**: transfer + native **full run merged** (4,957 annotated features)

- Transfer executed with default thresholds → 2,182/2,860 (76.3%) of EB's
  SIRIUS-annotated features assigned to local features; 1,773 local features
  annotated.
- Native SIRIUS pilot (150 targets → 149 usable spectra, 5 shards, job
  `27605104`) **completed and merged 2026-08-20**: 112 formula / 99 structure
  identifications folded in as an interim step (superseded by the full run
  below).
- **Native SIRIUS full run (job `27718540`, 127 shards, all 3,927 charge-1+
  targets) COMPLETED and MERGED 2026-08-25.** One shard (`_11`) hit the
  array's 1 h per-task time limit (`TIMEOUT`) and was resubmitted standalone
  on the `short` partition with `--time=02:00:00` (job `27748769_11`,
  completed in 1h07m) before merging. Merge + fold-in (`merged_full/` →
  `--native-label native-full-e9838293-bagel`) produced 3,184 formula /
  2,764 structure / 3,062 formula-CANOPUS / 2,764 structure-CANOPUS
  identifications, superseding the 112-row pilot (`native-EB97X`) for every
  feature it also covered. **Final: 4,957 local features annotated (1,773
  transferred + 3,184 native)**, 4,268 with a structure name. 0 new
  merge-conflicts from the fold-in.
- Downstream `differential-features-primary` / `feature-tables-primary`
  re-run 2026-08-25 with the full annotation set: 103,637 significant
  feature-rows, **2,639 bioactivity-flagged** (up from 1,106 with the pilot).
- Native pilot over-covered 6 features (`545`, `601`, `1905`, `2729`, `8042`,
  `8793`) — target selection did not exclude them, so each carried **both** a
  transferred row and a `native-EB97X` row before the full-run merge
  superseded the pilot rows for any feature the full run also covered.
- See "Known caveats" below for the merge-collapse conflicts.

## Datasets

- `EB/analysis/sirius_annotation/sirius_annotations.tsv` + `sirius_targets.csv`
  + `sirius_targets.mgf` — SIRIUS 6.3.12 outputs on the EB (MZmine3) feature
  table (read-only source). SIRIUS annotated 2,860/4,107 features (2,860
  formula, 2,647 structure); the other 886 are `2+`/`3+`/`4+` features that
  SIRIUS 6.3.12 formula tool skips ("Do not support multiple charges yet").
- `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/aligned_features_filled.mgf`
  — this project's feature spectra (`SCANS` == feature row id).
- `data/raw/gnps2_e9838293_bagel/nf_output/feature_finding/feature_finding_results/aligned_features.csv`
  — this project's feature table (component/row id, `feature_mz`, `feature_rt`).

## Method

### Key insight: how the transfer works

- **Join key**: EB `row ID` == `SCANS` in `sirius_targets.mgf`; this project's
  feature id == `SCANS` in `aligned_features_filled.mgf`.
- **Assign** each EB-annotated feature to a local feature via:
  1. exact text-match by precursor charge/adduct class (singletons matched
     directly by composition),
  2. `m/z` within `--ppm` **and** RT within `--rt` minutes, disambiguated by
     `--cos-min` MS2 spectral cosine (L2-normalized, fragment tolerance
     `--frag-tol` Da).

### Assignment decision tree (per EB feature)

1. Find local candidates: precursor `m/z` within `ppm`, RT within `rt`.
2. If a candidate's MS2 signature is empty (no peaks), it is kept as a
   pending candidate only if it is the sole composition match; otherwise
   filtered by MS2.
3. Scored (non-empty MS2) candidates get an MS2 cosine; pick best; require
   `cos >= cos_min`:
   - exactly one qualifying candidate → `ms2_winner`
   - `>1` tie on cos → `ms2_tie` (highest cosine), flagged for audit,
   - `0` candidates above threshold → `low_cosine`.
4. No candidates in the ppm/RT window → `no_match`.
5. No local candidate has MS2 at all → `no_ms2_candidate`.
6. If a candidate matches multiple EB features that **disagree on formula** →
   `merged_conflict` (see caveats).

### Defaults chosen

```
--ppm 10 --rt 0.5 --cos-min 0.7 --frag-tol 0.05
```

Threshold sweep (validated before running): (5 ppm, 0.3) → 2,044;
(15, 0.7) → 2,217; (20, 1.0) → 2,240. Defaults balance precision vs recovery.

## Results

| Metric | Value |
|--------|-------|
| EB SIRIUS-annotated features | 2,860 |
| Assigned (mapped to a local feature) | 2,182 (76.3%) |
| ─ `ms2_winner` | 472 |
| ─ `ms2_tie` | 533 |
| ─ unique formula/structure (single hits) | 1,177 |
| Not assigned | 678 (266 `no_match`, 50 `no_ms2_candidate`, 362 `low_cosine`) |
| **Local features annotated** (transferred only) | **1,773** |
| ─ with formula | 1,773 (100%) |
| ─ with structure name | 1,657 |

**After the full native run (2026-08-25):**

| Metric | Value |
|--------|-------|
| Native SIRIUS targets (charge 1+, un-annotated by transfer) | 3,927 |
| Usable native spectra (after degenerate-block filter) | 3,810 |
| Native formula identifications | 3,184 |
| Native structure identifications | 2,764 |
| **Local features annotated (final, accumulated)** | **4,957** |
| ─ transferred | 1,773 |
| ─ native | 3,184 |
| ─ with structure name | 4,268 |

## Outputs

| File | Description |
|------|-------------|
| `sirius_annotations.tsv` | Accumulated annotation table — **one row per local feature id**; the table to join onto `aligned_features.csv` (`row ID`). |
| `sirius_transfer_map.tsv` | Full transfer audit — one row per **EB source feature** (2,860 rows, incl. unassigned), so every assignment and every tie/conflict can be reconstructed. |
| `scripts/import_sirius_transfer.py` | Importer: runs the transfer, accumulates into `sirius_annotations.tsv`, and distills a future native SIRIUS merged dir into the same schema. |

### Column glossary (`sirius_annotations.tsv`)

- `row ID` — local feature id (join key to the feature table)
- `sirius_formula`, `sirius_adduct` — formula/adduct of the retained hit
- `sirius_structure_name`, `sirius_structure_smiles`, `sirius_structure_confidence`
- `sirius_npc_pathway`, `sirius_npc_class`, `sirius_classyfire_class` (CANOPUS)
- `sirius_source_feature_id`, `sirius_source_run` — which EB feature/run supplied it
- `annotation_origin` — `transferred` (from EB) or `native` (future run)
- `n_candidates`, `source_mz/source_rt`, `feature_mz/feature_rt`, `ppm_error`,
  `rt_delta_min`, `ms2_cosine`, `match_class`, `match_status`
- `n_sirius_hits` — number of distinct source feature ids collapsed onto this
  feature (concordant duplicate hits are expected)
- `sirius_hit_ids`, `sirius_hit_formulas` — semicolon-joined detail of all hits
- `merged_conflict` — **True** when the collapsed hits disagree on molecular
  formula (see caveats)

## Known caveats

1. **Merge-collapse conflicts (69 features).** The Everything-Bagel alignment
   is more permissive than the MZmine3 run, so several EB features can map to
   a *single* local feature — 296 features have >1 hit; of those, **69 have
   hits that disagree on molecular formula** (isobaric / co-isolated pairs,
   e.g. namalide d vs mycalamide A). The accumulated table keeps the
   highest-priority hit per feature; the discarded hit is fully visible in
   `sirius_transfer_map.tsv` (`sirius_hit_ids` / `sirius_hit_formulas`).
   Treat `merged_conflict == True` features as needing manual review.
2. MS2 cosine threshold 0.7 is permissive — always use `ms2_cosine` when
   ranking, not just the name.
3. EB used MZmine3 feature table; a compound absent there is absent here.

## Folding in a future native SIRIUS run

When SIRIUS is run directly on these features, run the merger to upgrade
`transferred` rows in place:

```bash
python3 scripts/import_sirius_transfer.py \
  --native-merged /path/to/sirius_results/merged \
  --native-label MYNATIVELABEL
```

- `--native-merged` must contain `formula_identifications.tsv` (required);
  `structure_identifications.tsv` and `canopus_*_summary.tsv` optional.
- Merging is **keyed by `row ID`** (native `mappingFeatureId` == local feature
  id). **Native rows win** over transferred rows per feature id
  (rank: `native` < `transferred`; then structure-hit; then
  `sirius_structure_confidence`).
- Transferred rows that the native run did not cover are kept.
- `--fresh` rebuilds `sirius_annotations.tsv` from scratch (ignores existing).
- If the native run overwrote annotations for a conflicting feature, the
  transferred conflict columns are recomputed after the merge.

(`argparse` is used, so `--help` lists every option, including the
`--ppm/--rt/--cos-min/--frag-tol` transfer thresholds.)

## Native SIRIUS run on the un-annotated features

The transfer left **3,927 local features** un-annotated but SIRIUS-runnable
(`has_ms2` and singly-charged). Native SIRIUS fills these in.

### Target set reduction

| Filter | Count |
|--------|-------|
| Features with MS2 | 6,453 |
| − already annotated (transferred) | −1,773 |
| **un-annotated with MS2** | 4,680 |
| − multiply charged (2+/3+; SIRIUS 6.3.12 refuses these) | −753 |
| **charge 1+ only (the native target set)** | **3,927** |

(Detection across analysis samples barely cuts this — MS2 features are all
detected broadly, median 59/90 use-in-analysis samples — so sample-presence is
not a useful reducer here. Charge state is the real lever.)

### Sharding / failure handling

The un-annotated set is split into small shards of ~20–40 features (default 30
per shard). Each shard is an **independent SIRIUS run into its own output
project space** (no shared project space), so a failed array task can be
re-submitted without touching the others. SIRIUS 6.3.12 serializes login-token
access, so the array runs strictly serial (`%1`) — the confirmed-safe
arrangement (EB / Rhodotorula experience). Per-feature runtime and shard size
can be tuned from the pilot results.

### Pipeline (scripts/run_sirius_native.sh)

```bash
# pilot: 150 targets (--max-features), 30 spectra/shard, 1 concurrent
bash analysis/sirius_annotation/scripts/run_sirius_native.sh 30 150 1
# full run: all charge-1+ targets
bash analysis/sirius_annotation/scripts/run_sirius_native.sh 30 0 1
```

The orchestrator runs, in order: **1)** `select_native_targets.py`
(`has_ms2` & `charge=1` & not in `sirius_annotations.tsv`; `--max-features`
draws a reproducible seed-sampled subset for pilots) → **2)**
`export_native_mgf.py` (pulls each target's MS2 block from
`aligned_features_filled.mgf`, keyed on `SCANS` == `row ID`) → **3)** shards the
MGF (EB's `shard_mgf.py`, round-robin, equal counts) → **4)** submits
`run_sirius_native.sbatch` to the `short` queue
(`--account=stajichlab`, 4 cpus, 16G, `SIRIUS_HEAP_GB=12`,
`formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 fingerprint canopus
structures write-summaries`) as a SLURM array, one shard per task.

**Degenerate-block guardrail**: `export_native_mgf.py` drops target blocks that
SIRIUS cannot use (missing/wrong `CHARGE` tag, `PEPMASS <= 0`, no positive
peaks) and reports them; e.g. the pilot dropped row ID `6813` (`CHARGE=0`,
`PEPMASS=0.0`) even though the feature table marks it `has_ms2`.

### Pilot (job 27605104)

150 targets (seed 1234) → 149 usable spectra → 5 shards of 30 (30/30/30/30/29).

**Completed 2026-08-20** (status above updated). Per-shard
`formula_identifications.tsv` rows: 24 / 11 / 24 / 24 / 29 → 112 merged
formula rows; 99 structure, 107 CANOPUS-formula, 99 CANOPUS-structure. All 5
array tasks finished comfortably within the `short`-queue limits (max walltime
< 24 h, heap 12G), validating shard size 30 as safe for the full run. No
denovo/spectral matches (0 rows) from the pilot's formula-first settings.

Merge + fold in (executed 2026-08-20):

```bash
python3 /bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/scripts/sirius_container_pipeline/merge_sirius_shards.py \
  --shard-root analysis/sirius_annotation/sirius_native_results \
  --out-dir analysis/sirius_annotation/sirius_native_results/merged
python3 analysis/sirius_annotation/scripts/import_sirius_transfer.py \
  --native-merged analysis/sirius_annotation/sirius_native_results/merged \
  --native-label native-EB97X
```

After fold-in: 1,885 assigned local features (1,773 transferred + 112 native),
0 new merge-conflicts.

### Full run (job 27718540)

3,927 charge-1+ targets → 3,810 usable spectra → 127 shards of 30
(`shard_000`–`shard_126`). Submitted via `bash
analysis/sirius_annotation/scripts/run_sirius_native.sh 30 0 1` to the
`short` partition, `%1` serial. 126/127 shards completed within the
per-task 1 h array time limit; shard `_11` hit `TIMEOUT` at 1h00m and was
resubmitted standalone with a longer time limit:

```bash
sbatch --partition=short --time=02:00:00 --array=11 \
  --export=ALL,SHARD_DIR=$PWD/analysis/sirius_annotation/shards_native,OUT_DIR=$PWD/analysis/sirius_annotation/sirius_native_results,SIF=/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif \
  analysis/sirius_annotation/scripts/run_sirius_native.sbatch
```

Completed 2026-08-25 (job `27748769_11`, 1h07m). Merge + fold in (executed
2026-08-25, after all 127 shards were COMPLETED):

```bash
python3 /bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/scripts/sirius_container_pipeline/merge_sirius_shards.py \
  --shard-root analysis/sirius_annotation/sirius_native_results \
  --out-dir analysis/sirius_annotation/sirius_native_results/merged_full
python3 analysis/sirius_annotation/scripts/import_sirius_transfer.py \
  --native-merged analysis/sirius_annotation/sirius_native_results/merged_full \
  --native-label native-full-e9838293-bagel
```

After fold-in: **4,957 assigned local features (1,773 transferred + 3,184
native)**, 0 new merge-conflicts. Downstream primary differential tables were
then regenerated (`pixi run differential-features-primary && pixi run
feature-tables-primary`): 103,637 significant feature-rows, 2,639
bioactivity-flagged.

**Note the output-dir naming**: the pilot's merge lives in `.../merged/`
(112 rows, `native-EB97X`, 2026-08-20) and the full run's merge lives in
`.../merged_full/` (3,184 rows, `native-full-e9838293-bagel`, 2026-08-25) —
do not reuse `merged/` for the full run's output, and do not confuse the two
labels when auditing `sirius_source_run` in `sirius_annotations.tsv`.

## Reproducibility

```bash
# transfer (defaults) + accumulate
python3 scripts/import_sirius_transfer.py
# merge a completed native run
python3 scripts/import_sirius_transfer.py --native-merged <dir> --native-label <label>
```

Full CLI help: `python3 scripts/import_sirius_transfer.py --help`.
Environment: repo pixi env
(`/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS/.pixi/envs/default/bin/python`).
Runtime is short (minutes); the transfer runs sequentially in a single process.
