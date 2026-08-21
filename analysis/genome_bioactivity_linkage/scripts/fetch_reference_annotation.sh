#!/usr/bin/bash -l
# analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh
#
# Downloads NCBI protein FASTA + GFF3 for the Bd JEL423 and Bsal AMFP13
# reference assemblies into analysis/genome_bioactivity_linkage/results/reference_annotation/.
# Idempotent: skips an assembly whose output already exists.
#
# Uses the raw NCBI datasets v2 REST download endpoint via curl (confirmed
# working: GET returns HTTP 200 with a valid zip, even though a `curl -I`
# HEAD request against the same URL returns 405 Method Not Allowed -- that
# 405 is expected and not a sign the endpoint is broken). No fallback to the
# ncbi-datasets-cli package was needed.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
mkdir -p "${OUT_ROOT}"

declare -A ASSEMBLIES=(
    [dendrobatidis]="GCA_048537975.1"
    [salamandrivorans]="GCA_002006685.2"
)

for species in "${!ASSEMBLIES[@]}"; do
    acc="${ASSEMBLIES[$species]}"
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/protein.faa" ]; then
        echo "SKIP ${species}: ${outdir}/protein.faa already exists"
        continue
    fi
    mkdir -p "${outdir}"
    zip="${outdir}/${acc}.zip"
    curl -sSL -o "${zip}" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=PROT_FASTA,GENOME_GFF"
    unzip -o -q "${zip}" -d "${outdir}"
    faa="$(find "${outdir}/ncbi_dataset/data" -iname 'protein.faa' | head -1)"
    gff="$(find "${outdir}/ncbi_dataset/data" -iname '*.gff' | head -1)"
    [ -n "${faa}" ] && cp "${faa}" "${outdir}/protein.faa"
    [ -n "${gff}" ] && cp "${gff}" "${outdir}/genomic.gff"
    if [ ! -s "${outdir}/protein.faa" ]; then
        echo "ERROR: no protein.faa downloaded for ${species} (${acc})" >&2
        exit 1
    fi
    echo "DONE ${species}: ${outdir}/protein.faa"
done
