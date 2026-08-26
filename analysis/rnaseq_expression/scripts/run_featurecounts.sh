#!/usr/bin/bash -l
# Count reads-per-gene on the completed STAR BAMs with featureCounts,
# using a gffread-converted GTF from the NCBI reference annotation
# (STAR --quantMode GeneCounts was useless because the raw NCBI GFF3 exon
# lines carry no gene_id attribute -- see RNASEQ_EXPRESSION.md).
# Runs all three strandedness modes; the highest-Assigned one is the
# correct library type for these public runs (unknown a priori).
#SBATCH -p short -N 1 -n 1 -c 8 --mem 24gb --time=02:00:00
#SBATCH --job-name=gbl_fc
#SBATCH --output=logs/gbl_featurecounts.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
source /etc/profile.d/modules.sh 2>/dev/null || true
module load subread/2.0.6

SPECIES_LIST=(dendrobatidis salamandrivorans)
species="${SPECIES_LIST[${SLURM_ARRAY_TASK_ID}]}"
SAMPLES="${REPO_DIR}/analysis/rnaseq_expression/scripts/samples.tsv"
gtf="${REPO_DIR}/analysis/rnaseq_expression/results/gtf/${species}.gtf"
outroot="${REPO_DIR}/analysis/rnaseq_expression/results/gene_counts/${species}"
mkdir -p "${outroot}" "${REPO_DIR}/logs"

# Derive BAM paths from samples.tsv (exact SRR ids) instead of a fragile glob --
# SRR ids are not a clean prefix family (e.g. SRR2768387{9,80,81},
# SRR13012113/117/121/125/129), and a sloppy pattern silently drops samples.
bams=$(awk -F '\t' -v s="${species}" '$2==s {print r"/analysis/rnaseq_expression/results/star_align/"$1"/Aligned.sortedByCoord.out.bam"}' r="${REPO_DIR}" "${SAMPLES}" | sort)
echo "${species}: $(echo ${bams} | wc -w) BAMs"

for s in 0 1 2; do
    out="${outroot}/counts_s${s}.txt"
    if [ -s "${out}.summary" ]; then
        echo "SKIP -s ${s}: ${out} exists"
        continue
    fi
    featureCounts -a "${gtf}" -o "${out}" -T "${SLURM_CPUS_PER_TASK:-8}" \
        -s "${s}" -p --countReadPairs ${bams} 2>"${outroot}/featurecounts_s${s}.log"
    echo "DONE -s ${s}: $(grep -P '^Assigned\t' ${out}.summary | cut -f2)"
done
echo "ALLDONE ${species}"
