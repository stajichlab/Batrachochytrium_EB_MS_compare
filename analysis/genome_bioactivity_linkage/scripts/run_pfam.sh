#!/usr/bin/bash -l
# Runs PFAM hmmsearch on both species' BFD-predicted proteomes locally,
# since BFD's own shared pfam_hmmscan run has no output at all for either
# Batrachochytrium genome as of 2026-08-25 (confirmed: zero files matching
# either locustag under BFD_ROOT/results/function/pfam_hmmscan). Mirrors
# BFD's own non-MPI hmmsearch invocation exactly (see
# BFD/Fungi_BFD_runs/nextflow/modules/BFD/PFAM/main.nf) so
# parse_pfam_domains.py needs no changes: same Pfam-A.hmm database,
# --cut_ga --noali, .domtblout/.tblout naming.
#
# Output: analysis/genome_bioactivity_linkage/results/pfam_hmmscan/<species>/<LOCUSTAG>.{domtblout,tblout}.gz
# (paths.find_bfd_output falls back here when BFD's own run has nothing).
#SBATCH -p short -N 1 -n 1 -c 8 --mem 16gb --time=02:00:00
#SBATCH --job-name=gbl_pfam
#SBATCH --output=logs/gbl_pfam.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/pfam_hmmscan"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load db-pfam
module load hmmer/3.4

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
domtbl="${outdir}/${tag}.domtblout"

if [ -s "${domtbl}.gz" ]; then
    echo "SKIP ${species}: ${domtbl}.gz already exists"
    exit 0
fi
fasta="${PROTEINS[$species]}"
[ -s "${fasta}" ] || { echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2; exit 1; }

mkdir -p "${outdir}"
hmmsearch --cut_ga --noali --cpu "${SLURM_CPUS_PER_TASK:-8}" \
    --domtbl "${domtbl}" \
    --tblout "${outdir}/${tag}.tblout" \
    "${PFAM_DB}/Pfam-A.hmm" "${fasta}" > /dev/null
gzip -f "${domtbl}" "${outdir}/${tag}.tblout"
echo "DONE ${species}: ${domtbl}.gz"
