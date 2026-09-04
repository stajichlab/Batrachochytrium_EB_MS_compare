# PEPTIDE_PROVENANCE — de novo fragment-tag matching

## Purpose

Settle H1 vs H2 for the blank-clearing, MS²-backed liquid-fraction peptides:

- **H1 (biosynthesis)** — products of the Tier-1 NRPS.
- **H2 (proteolysis)** — fragments of medium protein released by secreted
  fungal proteases (Bd 1% tryptone = whole-casein digest; Bsal 50% TGHL =
  tryptone + gelatin hydrolysate).

Two earlier attempts could not answer it. The proline-composition test is
**retracted** — it inherited the composition of SIRIUS's structure database.
Name-based sequence matching was negative but used database *guesses*, so a
true fragment named as the wrong isomer would be missed. This test uses the
**spectra**, so it is immune to database bias.

Repro: `pixi run fetch-substrates && pixi run fragment-tags`

## Method

1. Keep the `--top-peaks` most intense peaks; build a spectrum graph with an
   edge wherever two peaks differ by an amino-acid residue mass (±`--tol`).
2. Enumerate paths of `--tag-len` edges — each spells a de novo sequence tag.
3. Match tags as substrings of the six substrate proteins (αS1/αS2/β/κ casein,
   collagen α1(I)/α2(I); 3,677 aa), with **Leu≡Ile collapsed** (isobaric).
4. Score against a **letter-shuffled decoy tag set** (composition preserved).

### Three design choices forced by controls

Each is the *opposite* of the obvious approach, and each was caught by a
control rather than by reasoning:

**Match in both orientations.** A b/y spectrum spells the peptide twice — the
b-series N→C, the y-series C→N — so the tag set contains a sequence and its
reverse.

**Randomize the tag, not the substrate.** Because of the above, a
reversed-sequence decoy is structurally invalid: it scored exactly **1.00×**
against a known casein peptide at every tag length from 3 to 6. A
shuffled-*substrate* decoy is separately biased — collagen is a Gly-X-Y repeat
carrying only 0.43 distinct 3-mers per residue, so shuffling inflates the
decoy vocabulary by 20% and *penalises* the real sequence (first run: 0.86×).

**Use 5-residue tags, not 3.** Tag length sets how much of k-mer space the
substrate occupies:

| tag | distinct in substrate | possible (19 letters) | chance coverage |
|---|---|---|---|
| 3-mer | 1,468 | 6,859 | **21.4%** |
| 4-mer | 2,257 | 130,321 | 1.7% |
| 5-mer | 2,728 | 2,476,099 | **0.11%** |

At 3 residues a tag matches casein by luck one time in five, which is exactly
what the sweep shows (0.98–1.00× at every peak depth).

## Validation — the test has power

| control | 5-mer enrichment | p |
|---|---|---|
| **Positive** — synthetic b/y spectrum of β-casein 60–68 (`LQDKLHPFA`) + 50 noise peaks | **55.2×** | <0.001 |
| **Negative** — random 9-mer not in casein | **0.0×** | 1.000 |

Tag extraction itself is exact: 6 of 7 contiguous 3-mers of the control
peptide were recovered from its noisy spectrum. So a negative result below is
a real absence of signal, not a blind test.

## Result — NEGATIVE

Operating point `--top-peaks 150 --tag-len 5 --n-shuffle 500`:

| quantity | value |
|---|---|
| shortlist spectra yielding 5-mer tags | 46 / 267 (17%) |
| tag hits, real substrates | 72 |
| tag hits, composition-matched decoy | 46.1 |
| **aggregate enrichment** | **1.56×** |
| Wilcoxon (real > decoy across spectra) | **p = 0.15** |
| spectra at p < 0.05 | 4 / 46 (2.3 expected by chance), **binomial p = 0.20** |

**The shortlist is not detectably casein or collagen digest.** Against a
positive control that reaches 55×, an observed 1.56× that fails both
aggregate tests is an absence of signal. This independently confirms the
earlier name-based substring result through a completely different route.

Peak-depth / tag-length sweep (`sensitivity_sweep.tsv`,
`fragment_tag_sweep.png`) — the conclusion is stable across the grid:

| tag length | enrichment range across peak depths |
|---|---|
| 3 | 0.98 – 1.00× (no discrimination, as predicted) |
| 4 | 1.04 – 1.40× |
| 5 | 1.56 – 2.01× |

## The one directional hint

All four nominally-significant spectra are **Bsal, none Bd** — and only
Bsal's medium contains gelatin. The strongest, feature 943, is SIRIUS-named
`H-Pro-Leu-Glu-Pro-Ser-Gly-Gly-` — Pro/Gly-rich, i.e. collagen-like.

This is the same direction as the independent hydroxyproline-immonium result
(Hyp⁺ in 26.3% of Bsal liq peptide spectra vs 9.2% Bd, OR 3.53, p=7.1e-16).
But 4 significant of 46 against 2.3 expected is **not** significant on its own
(p=0.20), so it is a hint consistent with an independently-supported finding,
not evidence in its own right.

## What this means for the project goal

A minority of the Bsal liquid peptide pool is plausibly gelatin-derived. The
**bulk of the shortlist in both species is neither casein nor collagen
digest** — so the "secreted metabolome is medium proteolysis" reading, already
retracted on the composition side, is not rescued at the sequence level
either.

That leaves the shortlist's origin genuinely open: NRPS products, other
fungal peptides, or medium components too modified to match their parent
sequence. Distinguishing those still requires a defined or ¹³C/¹⁵N-labelled
medium — an experimental fix, not an analytical one.

## Caveats

- Only 17% of shortlist spectra yield 5-mer tags; the test is silent on the
  other 83%. A negative is therefore "no detectable digest signal in the
  fraction we can read", not "proven not digest".
- Tags are matched by sequence only, not anchored to a prefix mass —
  deliberately permissive, with the decoy null carrying the interpretation.
- Post-translational or process modifications (deamidation, oxidation) shift
  residue masses and would break a tag. Unmodified matching is the
  conservative choice and biases toward the negative observed here.
- A tag match is evidence about origin, not a structure assignment.

## Outputs

`fragment_tag_hits.tsv` (per spectrum), `sensitivity_sweep.tsv`,
`fragment_tag_sweep.png/.pdf`; substrates cached in
`reference_material/substrate_proteins/`.
