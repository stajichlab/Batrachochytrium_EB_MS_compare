#!/usr/bin/bash -l
# Runs MEROPS peptidase-family blastp on both species' BFD-predicted
# proteomes locally, since BFD's own shared merops run has no output at
# all for either Batrachochytrium genome as of 2026-08-26 (confirmed:
# zero files matching either locustag under
# BFD_ROOT/results/function/merops). Mirrors BFD's own invocation exactly
# (see BFD/Fungi_BFD_runs/nextflow/modules/BFD/MEROPS/main.nf) except
# using the host ncbi-blast/2.16.0+ module directly instead of BFD's
# containerized blastp -- same blast version, same flags, same
# merops_scan.lib database, so output format is identical.
#
# Output: analysis/genome_bioactivity_linkage/results/merops/<species>/<LOCUSTAG>.blasttab.gz
#SBATCH -p short -N 1 -n 1 -c 8 --mem 16gb --time=02:00:00
#SBATCH --job-name=gbl_merops
#SBATCH --output=logs/gbl_merops.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/merops"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load db-merops/124
module load ncbi-blast/2.16.0+

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
blasttab="${outdir}/${tag}.blasttab"

if [ -s "${blasttab}.gz" ]; then
    echo "SKIP ${species}: ${blasttab}.gz already exists"
    exit 0
fi
fasta="${PROTEINS[$species]}"
[ -s "${fasta}" ] || { echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2; exit 1; }

mkdir -p "${outdir}"
blastp -query "${fasta}" \
    -db "${MEROPS_DB}/merops_scan.lib" \
    -out "${blasttab}" \
    -num_threads "${SLURM_CPUS_PER_TASK:-8}" \
    -seg yes -soft_masking true \
    -max_target_seqs 10 \
    -evalue 1e-10 \
    -outfmt 6 \
    -use_sw_tback
gzip -f "${blasttab}"
echo "DONE ${species}: ${blasttab}.gz"
