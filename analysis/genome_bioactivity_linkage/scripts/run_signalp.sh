#!/usr/bin/bash -l
# Runs SignalP 6 on both species' BFD-predicted proteomes locally, since
# BFD's own shared signalp run has no output at all for either
# Batrachochytrium genome as of 2026-08-25. Mirrors BFD's own invocation
# (see BFD/Fungi_BFD_runs/nextflow/modules/BFD/SIGNALP/main.nf) -- flags
# match (-org euk --mode fast -format txt), so
# merge_secretion.load_signalp_gff3 needs no changes.
#
# Uses the GPU build (signalp/6.0h-gpu): the CPU build (job 27753328) ran
# at ~2.6 sequences/s and TIMEOUT'd on salamandrivorans (19,449 proteins,
# 92% done at the 2h limit) -- confirmed 2026-08-25.
#
# Output: analysis/genome_bioactivity_linkage/results/signalp/<species>/<LOCUSTAG>.signalp.gff3.gz
#SBATCH -p short_gpu -N 1 -n 1 -c 8 --mem 16gb --gres=gpu:1 --time=02:00:00
#SBATCH --job-name=gbl_signalp
#SBATCH --output=logs/gbl_signalp.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/signalp"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
# signalp/6 (default = 6.0i) is a broken install on this cluster missing
# the fast-mode model weights (distilled_model_signalp6.pt) -- confirmed
# 2026-08-25 (job 27753303 FAILED with FileNotFoundError). 6.0h has the
# complete model set (fast + sequential slow-mode models).
module load signalp/6.0h-gpu

declare -A PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)
declare -A LOCUSTAG=(
    [dendrobatidis]="FCC698BD"
    [salamandrivorans]="F61BA062"
)

SPECIES_LIST=(dendrobatidis salamandrivorans)
species="${SPECIES_LIST[${SLURM_ARRAY_TASK_ID}]}"
tag="${LOCUSTAG[$species]}"
outdir="${OUT_ROOT}/${species}"
gff3="${outdir}/${tag}.signalp.gff3"

if [ -s "${gff3}.gz" ]; then
    echo "SKIP ${species}: ${gff3}.gz already exists"
    exit 0
fi
fasta="${PROTEINS[$species]}"
[ -s "${fasta}" ] || { echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2; exit 1; }

mkdir -p "${outdir}"
workdir="${SCRATCH:?}/signalp_${species}"
rm -rf "${workdir}"
mkdir -p "${workdir}"
signalp6 -fasta "${fasta}" -od "${workdir}" -org euk --mode fast -format txt \
    --write_procs "${SLURM_CPUS_PER_TASK:-8}" -bs 100
cp "${workdir}/output.gff3" "${gff3}"
gzip -f "${gff3}"
cp "${workdir}/prediction_results.txt" "${outdir}/${tag}.signalp.results.txt"
gzip -f "${outdir}/${tag}.signalp.results.txt"
rm -rf "${workdir}"
echo "DONE ${species}: ${gff3}.gz"
