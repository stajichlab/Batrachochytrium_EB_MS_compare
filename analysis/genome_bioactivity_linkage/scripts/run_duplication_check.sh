#!/usr/bin/bash -l
# Duplication-rate check on the BFD-predicted proteomes (Bsal's 19,449-protein
# proteome vs Bd's 8,396 was flagged in GENOME_BIOACTIVITY_LINKAGE.md known
# caveat #3 as possibly inflated by over-prediction/gene-model duplication).
# Runs an all-vs-all DIAMOND self-blastp and flags, per query, any OTHER
# subject that is a near-duplicate (>= dup_pident identity AND >= dup_cov
# query coverage) -- the diagnostic the writeup recommends before trusting
# per-gene counts. Also reports the longest-protein distribution.
#SBATCH -p short -N 1 -n 1 -c 8 --mem 24gb --time=01:30:00
#SBATCH --job-name=gbl_dupcheck
#SBATCH --output=logs/gbl_dupcheck.%A_%a.log
#SBATCH --array=0-1
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-/bigdata/stajichlab/shared/projects/Chytrid/Bd_massspec/Batrachochytrium_MS}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/duplication_check"

declare -A PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)
SPECIES_LIST=(dendrobatidis salamandrivorans)
species="${SPECIES_LIST[${SLURM_ARRAY_TASK_ID}]}"
fasta="${PROTEINS[$species]}"
outdir="${OUT}/${species}"
mkdir -p "${outdir}" "${REPO_DIR}/logs"

source /etc/profile.d/modules.sh 2>/dev/null || true
module load diamond/2.1.7

[ -s "${fasta}" ] || { echo "ERROR: ${fasta} missing" >&2; exit 1; }

echo "${species}: $(grep -c '>' ${fasta}) proteins from ${fasta}"

# DIAMOND with --more-sensitive; self-vs-self, top 6 hits (self + up to 5 paralogs)
diamond makedb --in "${fasta}" -d "${outdir}/db" --quiet 2>"${outdir}/makedb.log" || true
if [ ! -s "${outdir}/db.dmnd" ]; then
    echo "makedb failed"; cat "${outdir}/makedb.log"; exit 1
fi
diamond blastp -d "${outdir}/db" -q "${fasta}" -o "${outdir}/self.blastp.tsv" \
    --outfmt 6 qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore \
    -k 6 -e 1e-5 --more-sensitive --threads "${SLURM_CPUS_PER_TASK:-8}" >/dev/null
echo "DONE ${species}: ${outdir}/self.blastp.tsv"

python3 - "${outdir}/self.blastp.tsv" <<'PY'
import sys
import pandas as pd
f = sys.argv[1]
cols = ["qseqid","sseqid","pident","length","qlen","slen","qstart","qend","sstart","send","evalue","bitscore"]
df = pd.read_csv(f, sep="\t", names=cols)
df = df[df["qseqid"] != df["sseqid"]]
df["query_cov"] = (df["qend"] - df["qstart"] + 1) / df["qlen"]
n_q = df["qseqid"].nunique()
# loose duplicate: >=90% identity over >=80% of the query
loose = df[(df["pident"] >= 90) & (df["query_cov"] >= 0.8)][["qseqid","sseqid","pident","query_cov","evalue"]].drop_duplicates("qseqid")
# strict duplicate: >=95% identity over >=90% query
strict = df[(df["pident"] >= 95) & (df["query_cov"] >= 0.9)][["qseqid","sseqid","pident","query_cov","evalue"]].drop_duplicates("qseqid")
# protein with the longest subject hit (paralog at >=50% pid) -- crude over-dup signal
print(f"query proteins with a non-self hit: {n_q}")
print(f"NEAR-DUPLICATES (>=90% id, >=80% query cov): {len(loose)} proteins have one")
print(f"STRICT DUPLICATES (>=95% id, >=90% query cov): {len(strict)} proteins have one")
if not loose.empty:
    top = loose.nlargest(15, "query_cov")
    print(top.to_string(index=False))
else:
    print("no near-duplicates found")
PY
echo "ALLDONE ${species}"
