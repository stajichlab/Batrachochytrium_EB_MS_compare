#!/usr/bin/bash -l
# Runs antiSMASH 8.0.4 (fungal mode) on the NCBI reference GBFF for both
# species, since BFD's own antiSMASH sub-run has not produced output for
# either genome (see GENOME_BIOACTIVITY_LINKAGE.md caveats). Uses the
# NCBI-annotated GBFF (fetch_reference_genbank.sh) rather than BFD's own
# gene models -- antiSMASH runs with --genefinding-tool none since the
# GBFF already carries CDS calls.
#
# IMPORTANT CAVEAT: this JSON's BGC-region record ids / coordinates are on
# NCBI's reference annotation, not BFD's predicted gene models. Downstream
# `parse_antismash_clusters.protein_in_bgc` compares BFD locus_tag genomic
# coordinates (from the BFD .gbk) against these regions -- that comparison
# is only valid if the NCBI assembly and the BFD-processed assembly share
# the same underlying sequence/contig coordinates. Verify contig/scaffold
# ids and lengths match between results/reference_annotation/<species>/genomic.gbff
# and BFD's .gbk (paths.bfd_gbk) before trusting has_bgc_context from this run.
#
# One SLURM array task per species (index 0=dendrobatidis, 1=salamandrivorans)
# rather than a single serial job -- on the `stajichlab` partition (30-day
# limit) since fungal --fullhmmer runtime on this cluster is not yet
# benchmarked and the `short` partition's 2h cap burned us once already on
# DeepTMHMM (see run_deeptmhmm.sh history).
#SBATCH -p stajichlab -N 1 -n 1 -c 8 --mem 16gb --time=08:00:00
#SBATCH --job-name=gbl_antismash_ref
#SBATCH --output=logs/gbl_antismash_ref.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
REF_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/antismash_ncbi"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load antismash/8.0.4

SPECIES_LIST=(dendrobatidis salamandrivorans)
declare -A OUT_NAME=(
    [dendrobatidis]="Batrachochytrium_dendrobatidis_JEL423"
    [salamandrivorans]="Batrachochytrium_salamandrivorans_AMFP13"
)

species="${SPECIES_LIST[${SLURM_ARRAY_TASK_ID}]}"
gbff="${REF_ROOT}/${species}/genomic.gbff"
outdir="${OUT_ROOT}/${species}"
json="${outdir}/${OUT_NAME[$species]}.json"

if [ -s "${json}" ]; then
    echo "SKIP ${species}: ${json} already exists"
    exit 0
fi
if [ ! -s "${gbff}" ]; then
    echo "ERROR: reference GBFF not found for ${species}: ${gbff}" >&2
    exit 1
fi
rm -rf "${outdir}"
mkdir -p "${outdir}"
antismash --taxon fungi --genefinding-tool none --fullhmmer \
    --cb-general --cb-knownclusters --cb-subclusters \
    -c "${SLURM_CPUS_PER_TASK:-8}" \
    --output-dir "${outdir}" \
    --output-basename "${OUT_NAME[$species]}" \
    "${gbff}"
echo "DONE ${species}: ${json}"
