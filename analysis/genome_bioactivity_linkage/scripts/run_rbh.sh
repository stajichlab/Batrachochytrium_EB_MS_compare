#!/usr/bin/bash -l
#SBATCH -p preempt -A preempt -N 1 -n 1 -c 8 --mem 16gb --time 0-02:00:00
#SBATCH --job-name=gbl_rbh
#SBATCH --output=logs/gbl_rbh.%j.log
#
# Reciprocal-best-hit DIAMOND blastp between each species' BFD gene models
# and its Task-4 reference annotation. Skips a species cleanly if its RBH
# output already exists.
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
REF_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/rbh"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load diamond

declare -A BFD_PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)

for species in "${!BFD_PROTEINS[@]}"; do
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/rbh.tsv" ]; then
        echo "SKIP ${species}: ${outdir}/rbh.tsv already exists"
        continue
    fi
    mkdir -p "${outdir}"
    bfd_fa="${BFD_PROTEINS[$species]}"
    ref_fa="${REF_ROOT}/${species}/protein.faa"
    if [ ! -s "${ref_fa}" ]; then
        echo "ERROR: reference protein.faa not found for ${species}: ${ref_fa} (run fetch_reference_annotation.sh first)" >&2
        exit 1
    fi
    diamond makedb --in "${bfd_fa}" -d "${outdir}/bfd_db"
    diamond makedb --in "${ref_fa}" -d "${outdir}/ref_db"
    diamond blastp -q "${bfd_fa}" -d "${outdir}/ref_db" -o "${outdir}/fwd.tsv" \
        --threads 8 --max-target-seqs 5 --evalue 1e-10
    diamond blastp -q "${ref_fa}" -d "${outdir}/bfd_db" -o "${outdir}/rev.tsv" \
        --threads 8 --max-target-seqs 5 --evalue 1e-10
    echo "DONE ${species}: ${outdir}/fwd.tsv, ${outdir}/rev.tsv"
done
