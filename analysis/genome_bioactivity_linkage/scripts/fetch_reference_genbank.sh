#!/usr/bin/bash -l
# analysis/genome_bioactivity_linkage/scripts/fetch_reference_genbank.sh
#
# Downloads NCBI's full annotated GenBank flat file (GENOME_GBFF: sequence +
# gene/CDS features in one file) for the Bd JEL423 and Bsal AMFP13 reference
# assemblies -- antiSMASH's preferred input format, since it carries gene
# calls directly (no separate genefinding step needed).
#
# This is a companion to fetch_reference_annotation.sh (which only pulls
# protein.faa + genomic.gff for the RBH cross-reference step); the GBFF
# fetched here feeds run_antismash_reference.sh instead. Idempotent:
# skips an assembly whose GBFF already exists.
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
    gbff="${outdir}/genomic.gbff"
    if [ -s "${gbff}" ]; then
        echo "SKIP ${species}: ${gbff} already exists"
        continue
    fi
    mkdir -p "${outdir}"
    zip="${outdir}/${acc}_gbff.zip"
    curl -sSL -o "${zip}" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=GENOME_GBFF"
    unzip -o -q "${zip}" -d "${outdir}"
    found="$(find "${outdir}/ncbi_dataset/data" -iname '*.gbff' | head -1)"
    if [ -z "${found}" ]; then
        echo "ERROR: no .gbff downloaded for ${species} (${acc})" >&2
        exit 1
    fi
    cp "${found}" "${gbff}"
    rm -f "${zip}"
    echo "DONE ${species}: ${gbff}"
done
