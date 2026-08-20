#!/usr/bin/env bash
# SIRIUS native-annotation chain as a SLURM array job on the short queue.
# One shard per array task. Concurrency is controlled by the --array=..%N flag
# on the sbatch command line in run_sirius_native.sh (default %1, strict serial
# -- SIRIUS 6.3.12 serializes login-token access and parallel arrays broke
# login in the sibling Rhodotorula project; %1 is the confirmed-safe default).
#
# Not submitted directly; run_sirius_native.sh sets SHARD_DIR/OUT_DIR/SIF and
# passes --array=0-N%1 on the sbatch command line.
#SBATCH --job-name=sirius_native_array
#SBATCH --account=stajichlab
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%A_%a.log
#SBATCH --error=%x-%A_%a.log
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"

: "${SHARD_DIR:?set via --export}"
: "${OUT_DIR:?set via --export}"
: "${SIF:?set via --export}"
: "${SLURM_ARRAY_TASK_ID:?set by SLURM array}"

module load singularity

SIRIUS_HEAP_GB="${SIRIUS_HEAP_GB:-12}"

SHARD_ID=$(printf 'shard_%03d' "${SLURM_ARRAY_TASK_ID}")
SHARD="${SHARD_DIR}/${SHARD_ID}.mgf"
SHARD_OUT="${OUT_DIR}/${SHARD_ID}"

echo "=== ${SHARD_ID} ==="

if [[ ! -f "$SHARD" ]]; then
  echo "Shard $SHARD not found -- nothing to do"
  exit 0
fi

mkdir -p "$OUT_DIR"

JAVA_OPTS="-Xmx${SIRIUS_HEAP_GB}G -XX:AOTMode=off" singularity exec --bind /bigdata,/scratch "$SIF" \
  sirius --cores "${SLURM_CPUS_PER_TASK:-4}" \
  --input "$SHARD" --output "$SHARD_OUT" \
  formula --ppm-max 15 --ppm-max-ms2 15 --candidates 10 \
  fingerprint \
  canopus \
  structures \
  write-summaries

echo "${SHARD_ID} completed"
