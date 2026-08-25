#!/usr/bin/bash -l
# Aligns each RNA-seq run (samples.tsv) against its species' STAR index,
# producing a coordinate-sorted BAM + per-gene counts (--quantMode
# GeneCounts, all 3 strandedness columns -- library strandedness for these
# public runs is not known a priori, so all 3 are kept and the correct
# column is picked during analysis, not here).
#
# Submit with --dependency=afterok:<download_jobid>:<index_jobid> so this
# never races an incomplete index or fastq download.
#SBATCH -p short -N 1 -n 1 -c 8 --mem 16gb --time=02:00:00
#SBATCH --job-name=gbl_rnaseq_align
#SBATCH --output=logs/gbl_rnaseq_align.%A_%a.log
#SBATCH --array=0-7
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
SAMPLES="${REPO_DIR}/analysis/rnaseq_expression/scripts/samples.tsv"
FASTQ_ROOT="${REPO_DIR}/data/raw/rnaseq_srp291769"
INDEX_ROOT="${REPO_DIR}/analysis/rnaseq_expression/results/star_index"
OUT_ROOT="${REPO_DIR}/analysis/rnaseq_expression/results/star_align"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load star/2.7.11b
module load samtools/1.19.2

line_num=$(( SLURM_ARRAY_TASK_ID + 2 ))
read -r srr species rep < <(sed -n "${line_num}p" "${SAMPLES}")
[ -n "${srr}" ] || { echo "ERROR: no sample at line ${line_num} of ${SAMPLES}" >&2; exit 1; }

fq_dir="${FASTQ_ROOT}/${srr}"
index_dir="${INDEX_ROOT}/${species}"
outdir="${OUT_ROOT}/${srr}"
bam="${outdir}/Aligned.sortedByCoord.out.bam"

if [ -s "${bam}" ]; then
    echo "SKIP ${srr} (${species} rep${rep}): ${bam} already exists"
    exit 0
fi
[ -s "${index_dir}/SAindex" ] || { echo "ERROR: STAR index missing for ${species}: ${index_dir}" >&2; exit 1; }

r1="${fq_dir}/${srr}_1.fastq.gz"
r2="${fq_dir}/${srr}_2.fastq.gz"
if [ -s "${r1}" ] && [ -s "${r2}" ]; then
    reads=("${r1}" "${r2}")
elif [ -s "${fq_dir}/${srr}.fastq.gz" ]; then
    reads=("${fq_dir}/${srr}.fastq.gz")
else
    echo "ERROR: no fastq.gz found for ${srr} in ${fq_dir}" >&2
    exit 1
fi

rm -rf "${outdir}"
mkdir -p "${outdir}"
STAR --runMode alignReads \
    --genomeDir "${index_dir}" \
    --readFilesIn "${reads[@]}" \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --quantMode GeneCounts \
    --outFileNamePrefix "${outdir}/" \
    --runThreadN "${SLURM_CPUS_PER_TASK:-8}"
samtools index "${bam}"
echo "DONE ${srr} (${species} rep${rep}): ${bam}"
