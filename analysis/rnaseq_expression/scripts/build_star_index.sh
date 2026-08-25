#!/usr/bin/bash -l
# Builds a STAR genome index per species from the NCBI reference genome
# (genome FASTA extracted from the already-downloaded genomic.gbff;
# annotation from the already-downloaded genomic.gff -- both fetched by
# fetch_reference_annotation.sh / fetch_reference_genbank.sh for the
# genome_bioactivity_linkage pipeline, reused here so the RNA-seq gene ids
# match the same NCBI annotation used for the RBH cross-reference step).
#
# Small fungal genomes need --genomeSAindexNbases reduced from STAR's
# default (14) per STAR's own manual formula: min(14, log2(genomeLength)/2 - 1).
# Computed per species below rather than hardcoded.
#SBATCH -p short -N 1 -n 1 -c 8 --mem 32gb --time=01:00:00
#SBATCH --job-name=gbl_star_index
#SBATCH --output=logs/gbl_star_index.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
REF_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
OUT_ROOT="${REPO_DIR}/analysis/rnaseq_expression/results/star_index"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load star/2.7.11b

SPECIES_LIST=(dendrobatidis salamandrivorans)
species="${SPECIES_LIST[${SLURM_ARRAY_TASK_ID}]}"

gbff="${REF_ROOT}/${species}/genomic.gbff"
gff="${REF_ROOT}/${species}/genomic.gff"
fasta="${REF_ROOT}/${species}/genomic.fna"
outdir="${OUT_ROOT}/${species}"

if [ -s "${outdir}/SAindex" ]; then
    echo "SKIP ${species}: ${outdir}/SAindex already exists"
    exit 0
fi
[ -s "${gbff}" ] || { echo "ERROR: missing ${gbff}" >&2; exit 1; }
[ -s "${gff}" ] || { echo "ERROR: missing ${gff}" >&2; exit 1; }

if [ ! -s "${fasta}" ]; then
    python3 - "${gbff}" "${fasta}" <<'PYEOF'
import sys
from Bio import SeqIO
gbff, fasta = sys.argv[1], sys.argv[2]
n = SeqIO.convert(gbff, "genbank", fasta, "fasta")
print(f"extracted {n} sequences to {fasta}")
PYEOF
fi

genome_len=$(python3 -c "
from Bio import SeqIO
print(sum(len(r) for r in SeqIO.parse('${fasta}', 'fasta')))
")
sa_index_nbases=$(python3 -c "
import math
n = min(14, int(math.log2(${genome_len}) / 2 - 1))
print(n)
")
echo "${species}: genome length ${genome_len} bp -> --genomeSAindexNbases ${sa_index_nbases}"

rm -rf "${outdir}"
mkdir -p "${outdir}"
STAR --runMode genomeGenerate \
    --genomeDir "${outdir}" \
    --genomeFastaFiles "${fasta}" \
    --sjdbGTFfile "${gff}" \
    --sjdbGTFtagExonParentTranscript Parent \
    --sjdbGTFfeatureExon exon \
    --sjdbOverhang 100 \
    --genomeSAindexNbases "${sa_index_nbases}" \
    --runThreadN "${SLURM_CPUS_PER_TASK:-8}"
echo "DONE ${species}: ${outdir}"
