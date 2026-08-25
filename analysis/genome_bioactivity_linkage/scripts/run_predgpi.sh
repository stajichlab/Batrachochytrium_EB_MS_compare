#!/usr/bin/bash -l
# Runs PredGPI on both species' BFD-predicted proteomes locally, since
# BFD's own shared predgpi run has no output at all for either
# Batrachochytrium genome as of 2026-08-25. Mirrors BFD's own invocation
# exactly (see BFD/Fungi_BFD_runs/nextflow/modules/BFD/PREDGPI/main.nf),
# so merge_secretion.load_predgpi_gff3 needs no changes.
#
# Output: analysis/genome_bioactivity_linkage/results/predgpi/<species>/<LOCUSTAG>.predgpi.gff3.gz
#SBATCH -p short -N 1 -n 1 -c 4 --mem 8gb --time=02:00:00
#SBATCH --job-name=gbl_predgpi
#SBATCH --output=logs/gbl_predgpi.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/predgpi"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load predgpi

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
gff3="${outdir}/${tag}.predgpi.gff3"

if [ -s "${gff3}.gz" ]; then
    echo "SKIP ${species}: ${gff3}.gz already exists"
    exit 0
fi
fasta="${PROTEINS[$species]}"
[ -s "${fasta}" ] || { echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2; exit 1; }

mkdir -p "${outdir}"
predgpi.py -f "${fasta}" -m gff3 -o "${gff3}"
gzip -f "${gff3}"
echo "DONE ${species}: ${gff3}.gz"
