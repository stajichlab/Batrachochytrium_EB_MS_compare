#!/usr/bin/env bash
# Orchestrator: build the native SIRIUS target set, export MGF, shard, and
# submit the equally-sized shards to the short queue as a SLURM array.
#
# This is the pilot-then-scale driver:
#   * default = PILOT: 150 features (charge 1+ only), 30 spectra/shard
#     -> ~5 shards, benchmarkable, restartable per shard
#   * full run  = --max-features 0  (0 = all 3,927 charge-1+ targets)
#
# Outputs (all under analysis/sirius_annotation/):
#   sirius_native_targets.csv        target list (row ID / m/z / RT / provenance)
#   sirius_native_targets.mgf        per-feature MS2 MGF pulled from the feature MGF
#   shards_native/shard_%03d.mgf     sharded MGFs (one per array task)
#   sirius_native_results/shard_%03d/  per-shard SIRIUS project space
#
# Each shard is an independent SIRIUS run into its own output project space, so
# a failed array task can be re-run/resubmitted without touching the others --
# that is the failure/restart model you asked for (no shared project space).
#
# Merge + import (after the array completes):
#   python3 /bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/scripts/sirius_container_pipeline/merge_sirius_shards.py \
#       --shard-root "$PWD/analysis/sirius_annotation/sirius_native_results" \
#       --out-dir "$PWD/analysis/sirius_annotation/sirius_native_results/merged"
#   python3 scripts/import_sirius_transfer.py \
#       --native-merged analysis/sirius_annotation/sirius_native_results/merged \
#       --native-label native-EB97X
#
# Usage:
#   scripts/run_sirius_native.sh [spectra_per_shard] [max_features] [concurrency]
#       (defaults: 30 150 1)
set -euo pipefail
# not a git repo; anchor to this file's location (repo root = 3 levels up)
cd "$(dirname "$0")/../../.."

SPECTRA_PER_SHARD="${1:-30}"
MAX_FEATURES="${2:-150}"                                  # 0 => all targets
CONCURRENCY="${3:-1}"                                     # %N parallel tasks
SIF=/bigdata/stajichlab/shared/singularity/sirius-6.3.12-linux-x64.sif
PY=python3
PIPELINE=analysis/sirius_annotation/scripts
MERGE=/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/scripts/sirius_container_pipeline/merge_sirius_shards.py

[[ -f "$SIF" ]] || { echo "Container image not found: $SIF" >&2; exit 1; }

SELECT_FLAGS=()
[[ "$MAX_FEATURES" =~ ^[0-9]+$ ]] && [[ "$MAX_FEATURES" -gt 0 ]] && SELECT_FLAGS=(--max-features "$MAX_FEATURES")

# Step 1: select targets (has_ms2 & charge 1+ & un-annotated)
$PY "$PIPELINE/select_native_targets.py" "${SELECT_FLAGS[@]:-}"

# Step 2: export MGF for the targets from the Everything-Bagel feature MGF
$PY "$PIPELINE/export_native_mgf.py" --charge 1

# Step 3: shard the MGF (round-robin, equal counts -- like EB's shard_mgf.py)
MGF=analysis/sirius_annotation/sirius_native_targets.mgf
SHARD_DIR=analysis/sirius_annotation/shards_native
N_SPECTRA=$(grep -c '^BEGIN IONS' "$MGF")
N_SHARDS=$(( (N_SPECTRA + SPECTRA_PER_SHARD - 1) / SPECTRA_PER_SHARD ))
echo "Sharding $N_SPECTRA spectra into $N_SHARDS shards (~${SPECTRA_PER_SHARD}/shard)"
rm -rf "$SHARD_DIR"
$PY /bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/EB/scripts/sirius_container_pipeline/shard_mgf.py \
  --input "$MGF" --out-dir "$SHARD_DIR" --n-shards "$N_SHARDS"

# Step 4: submit array job (one shard per task, %CONCURRENCY at a time)
OUT_DIR=analysis/sirius_annotation/sirius_native_results
mkdir -p analysis/sirius_annotation/logs

JOBID=$(sbatch --chdir="$PWD" \
  --array=0-$((N_SHARDS-1))%"$CONCURRENCY" \
  --export=ALL,SHARD_DIR="$PWD/$SHARD_DIR",OUT_DIR="$PWD/$OUT_DIR",SIF="$SIF" \
  --output="analysis/sirius_annotation/logs/%A_%a.log" \
  --error="analysis/sirius_annotation/logs/%A_%a.log" \
  --parsable "$PIPELINE/run_sirius_native.sbatch")

echo "Submitted array job $JOBID ($N_SHARDS shards, $CONCURRENCY at a time, short queue)"
echo "After it completes, merge and import with:"
echo "  $PY $MERGE --shard-root $PWD/$OUT_DIR --out-dir $PWD/$OUT_DIR/merged"
echo "  $PY $PIPELINE/import_sirius_transfer.py --native-merged $PWD/$OUT_DIR/merged --native-label native-EB97X"
