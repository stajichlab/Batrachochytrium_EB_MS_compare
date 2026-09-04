#!/usr/bin/env bash
# Fetch the growth-medium substrate proteins from UniProt.
#
# Bd medium  = 1% tryptone            -> pancreatic digest of WHOLE casein
#                                        (alpha-S1, alpha-S2, beta, kappa)
# Bsal medium = 50% TGHL              -> tryptone + GELATIN hydrolysate
#                                        (denatured collagen type I)
#
# Whole casein is the right reference, not beta-casein alone -- using the one
# fraction that happened to fit was an error in the retracted proline test.
set -euo pipefail
OUT="$(cd "$(dirname "$0")/../../.." && pwd)/reference_material/substrate_proteins"
mkdir -p "$OUT"
declare -A ACC=(
  [P02662]="casein_alphaS1" [P02663]="casein_alphaS2"
  [P02666]="casein_beta"    [P02668]="casein_kappa"
  [P02453]="collagen_a1_I"  [P02465]="collagen_a2_I"
)
for acc in "${!ACC[@]}"; do
  f="$OUT/${ACC[$acc]}.${acc}.fasta"
  if [ -s "$f" ]; then echo "  skip (present): $(basename "$f")"; continue; fi
  curl -sS --fail --retry 3 "https://rest.uniprot.org/uniprotkb/${acc}.fasta" -o "$f"
  n=$(grep -v '^>' "$f" | tr -d '\n' | wc -c)
  echo "  fetched ${acc} ${ACC[$acc]} (${n} aa)"
done
echo "substrates in $OUT"
