#!/usr/bin/bash -l
# Downloads the 8 RNA-seq runs listed in samples.tsv (3 Bd JEL423 reps +
# 5 Bsal AMFP13 reps, SRP291769) via sratoolkit prefetch + fasterq-dump.
# One SLURM array task per run. Confirmed compute nodes on the `short`
# partition have outbound internet access to NCBI (2026-08-25 nettest).
#
# NOTE: these RNA-seq runs are NOT from the liquid-culture growth
# condition sampled by this project's mass-spec data -- see the caveat in
# RNASEQ_EXPRESSION.md. They provide presence/absence expression evidence
# only, not condition-matched co-expression evidence.
#SBATCH -p short -N 1 -n 1 -c 4 --mem 8gb --time=02:00:00
#SBATCH --job-name=gbl_rnaseq_dl
#SBATCH --output=logs/gbl_rnaseq_dl.%A_%a.log
#SBATCH --array=0-7
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
SAMPLES="${REPO_DIR}/analysis/rnaseq_expression/scripts/samples.tsv"
OUT_ROOT="${REPO_DIR}/data/raw/rnaseq_srp291769"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load sratoolkit/3.2.0

# samples.tsv has a header; array index 0 -> data line 1 (tail -n +2, then +1 for 1-based sed)
line_num=$(( SLURM_ARRAY_TASK_ID + 2 ))
read -r srr species rep < <(sed -n "${line_num}p" "${SAMPLES}")
if [ -z "${srr}" ]; then
    echo "ERROR: no sample at line ${line_num} of ${SAMPLES}" >&2
    exit 1
fi

outdir="${OUT_ROOT}/${srr}"
if compgen -G "${outdir}"'/*.fastq.gz' > /dev/null; then
    echo "SKIP ${srr} (${species} rep${rep}): fastq.gz already present in ${outdir}"
    exit 0
fi
mkdir -p "${outdir}"

echo "=== ${srr} (${species} rep${rep}) ==="
prefetch --max-size 100G -O "${SCRATCH:?}" "${srr}"
fasterq-dump --split-files --threads "${SLURM_CPUS_PER_TASK:-4}" \
    -O "${outdir}" "${SCRATCH:?}/${srr}/${srr}.sra"
gzip -f "${outdir}"/*.fastq
rm -rf "${SCRATCH:?}/${srr}"
echo "DONE ${srr}: $(ls "${outdir}"/*.fastq.gz | wc -l) fastq.gz file(s)"
