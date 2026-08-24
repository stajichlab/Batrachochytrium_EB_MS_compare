#!/usr/bin/bash -l
# -p short_gpu (rather than preempt_gpu) is a DELIBERATE override of the
# partition, per the project owner's explicit instruction for this job --
# do not "fix" it back to preempt_gpu. -A preempt is still the project
# owner's stated blanket account preference for all SLURM jobs.
#SBATCH -p gpu -N 1 -n 1 -c 4 --mem 16gb --gres=gpu:1
#SBATCH --job-name=gbl_deeptmhmm
#SBATCH --output=logs/gbl_deeptmhmm.%j.log
#
# Runs DeepTMHMM on both Bd JEL423 and Bsal AMFP13 proteomes, writing
# TMRs.gff3 into analysis/genome_bioactivity_linkage/results/deeptmhmm/<species>/.
# Skips a species cleanly if its output already exists (Task 2's caching
# requirement from the spec).
#
# set -e: a real DeepTMHMM crash (predict.py raising, e.g. the CRF-decode
# ValueError seen on Bsal's one 4,777-aa outlier protein -- see
# EXCLUDED_PROTEIN_IDS below) must fail this job loudly, not print a false
# "DONE" and exit 0. Confirmed necessary: run 27718539 (2026-08-23) crashed
# mid-Bsal but the unset-`-e` version of this script still logged "DONE
# salamandrivorans" and the SLURM job reported COMPLETED, even though no
# TMRs.gff3 was ever written.
set -euo pipefail

SIF="${DEEPTMHMM_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif}"
# BASH_SOURCE resolves to empty/wrong paths in SLURM work directories --
# never use it here (see ~/.claude/CLAUDE.md). REPO_DIR must come from
# SLURM_SUBMIT_DIR (set by sbatch) or this hardcoded absolute fallback.
REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/deeptmhmm"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer

declare -A PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)

# Proteins excluded from DeepTMHMM input (space-separated IDs per species).
# F61BA062_016014-T1 (Bsal AMFP13, 4,777 aa -- the single longest protein in
# either proteome) crashes DeepTMHMM's topology-decoding step with
# "ValueError: the first two dimensions of emissions and tags must match,
# got (4777, 1) and (1, 1)" (run 27718539, 2026-08-23 log). Excluded here so
# the rest of the proteome can be processed; the excluded record is saved to
# ${OUT_ROOT}/<species>_excluded_proteins.fasta for separate investigation
# (plausibly a fused/misassembled gene model -- see the Bsal protein-count
# anomaly caveat in GENOME_BIOACTIVITY_LINKAGE.md).
declare -A EXCLUDED_PROTEIN_IDS=(
    [salamandrivorans]="F61BA062_016014-T1"
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

    excluded_ids="${EXCLUDED_PROTEIN_IDS[$species]:-}"
    if [ -n "${excluded_ids}" ]; then
        excluded_fasta="${OUT_ROOT}/${species}_excluded_proteins.fasta"
        filtered_fasta="${SCRATCH:?}/${species}.deeptmhmm_input.proteins.fa"
        awk -v ids="${excluded_ids}" -v keepf="${filtered_fasta}" -v exclf="${excluded_fasta}" '
            BEGIN { n = split(ids, arr, " "); for (i = 1; i <= n; i++) idset[">" arr[i]] = 1 }
            /^>/ { out = ($1 in idset) ? exclf : keepf }
            { print > out }
        ' "${fasta}"
        echo "EXCLUDED from ${species} DeepTMHMM input: ${excluded_ids} (saved to ${excluded_fasta})"
        fasta="${filtered_fasta}"
    fi

    rm -rf "${outdir}"
    apptainer exec --nv -B "${REPO_DIR}" -B "${BFD_ROOT}" -B "${SCRATCH:-/tmp}" "${SIF}" \
        bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta ${fasta} --output-dir ${outdir}"
    echo "DONE ${species}: ${outdir}/TMRs.gff3"
done
