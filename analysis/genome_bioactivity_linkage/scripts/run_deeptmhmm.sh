#!/usr/bin/bash -l
#SBATCH -p short_gpu -N 1 -n 1 -c 4 --mem 16gb --gres=gpu:1 --time 0-02:00:00
#SBATCH --job-name=gbl_deeptmhmm
#SBATCH --output=logs/gbl_deeptmhmm.%j.log
#
# Runs DeepTMHMM on both Bd JEL423 and Bsal AMFP13 proteomes, writing
# TMRs.gff3 into analysis/genome_bioactivity_linkage/results/deeptmhmm/<species>/.
# Skips a species cleanly if its output already exists (Task 2's caching
# requirement from the spec).
set -euo pipefail

SIF="${DEEPTMHMM_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif}"
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/deeptmhmm"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer

declare -A PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)

for species in "${!PROTEINS[@]}"; do
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/TMRs.gff3" ]; then
        echo "SKIP ${species}: ${outdir}/TMRs.gff3 already exists"
        continue
    fi
    fasta="${PROTEINS[$species]}"
    if [ ! -s "${fasta}" ]; then
        echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2
        exit 1
    fi
    rm -rf "${outdir}"
    apptainer exec --nv -B "${REPO_DIR}" -B "${BFD_ROOT}" "${SIF}" \
        bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta ${fasta} --output-dir ${outdir}"
    echo "DONE ${species}: ${outdir}/TMRs.gff3"
done
