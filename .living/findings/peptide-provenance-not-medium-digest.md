## F-006: The blank-clearing liquid peptides are NOT detectably casein/collagen digest — de novo fragment tags say so with validated power
**Status:** supported (negative result, positive-control-validated)
**Claim:** De novo MS2 sequence-tag matching of the 267 blank-clearing, MS2-backed liq shortlist spectra against the six growth-medium substrate proteins (bovine alphaS1/alphaS2/beta/kappa casein + collagen alpha1(I)/alpha2(I), 3,677 aa) gives **aggregate enrichment 1.56x over a composition-matched decoy, Wilcoxon p=0.15, and 4/46 spectra at p<0.05 against 2.3 expected by chance (binomial p=0.20)** — i.e. no signal. The test is not blind: a synthetic b/y spectrum of beta-casein 60-68 (`LQDKLHPFA`) plus 50 noise peaks scores **55.2x (p<0.001)** and a random 9-mer scores **0.0x (p=1.000)**, and tag extraction recovers 6 of 7 contiguous 3-mers of the control peptide. The conclusion is stable across a peak-depth x tag-length sweep (3-mers 0.98-1.00x, 4-mers 1.04-1.40x, 5-mers 1.56-2.01x). This independently confirms, from spectra rather than database names, the earlier negative name-based substring test.
**Implications:** The "secreted metabolome is dominated by proteolysis of medium protein" reading — already retracted on the proline-composition side — is **not rescued at the sequence level either**. The bulk of the shortlist in both species is neither casein nor collagen digest, leaving its origin genuinely open (NRPS products, other fungal peptides, or medium components too modified to match their parent sequence). One directional hint survives: all four nominally-significant spectra are Bsal and none Bd, and only Bsal's medium (50% TGHL) contains gelatin; the strongest (feature 943) is SIRIUS-named `H-Pro-Leu-Glu-Pro-Ser-Gly-Gly-`, Pro/Gly-rich and collagen-like. That is the same direction as the independent hydroxyproline-immonium contrast (Bsal 26.3% vs Bd 9.2%, OR 3.53, p=7.1e-16) but is not significant alone. Resolving the shortlist's origin requires a defined or 13C/15N-labelled medium — experimental, not analytical.
**Tags:** metabolomics, peptide-provenance, fragment-tags, de-novo-sequencing, negative-result, positive-control, casein, collagen, secreted-compounds, everything-bagel

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-03 | adduct-fix-and-provenance-2026-09-03 | 267 shortlist MS2 spectra + 6 UniProt substrates | Batrachochytrium_EB_MS_compare | aggregate 1.56x, Wilcoxon p=0.15, 4/46 spectra p<0.05 vs 2.3 expected (binom p=0.20); positive control 55.2x, negative control 0.0x | refutes H2 for the bulk of the shortlist |

### Open Questions
- Would allowing common modifications (deamidation, oxidation) on tag masses recover a digest signal that unmodified matching misses?
- The 83% of shortlist spectra that yield no 5-mer tag are untested — are they too sparse, chimeric, or non-peptidic?
- Do the 4 Bsal-only hits survive a targeted, prefix-mass-anchored re-analysis against collagen specifically?

---
