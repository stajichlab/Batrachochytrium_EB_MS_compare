# MOLECULAR_NETWORK

## Purpose
Molecular-family (component) tier for the FBMN network — GOALS.md goals 1
(characterize molecular families) and 4 (probe DeltaMZ ladders), which had
never been attempted before 2026-09-02.

This tier exists because per-feature evidence in this dataset is weak by
construction: n=5 per group, a BH threshold that only fires when ~16% of
features separate simultaneously, and a medium background that mimics the
biological hypothesis. A molecular family is a stronger unit of evidence —
a homologous series that is *entirely* blank-clearing is much harder to
attribute to medium background than any single feature.

## Method
1. Connected components from `filtered_pairs.tsv` (MS² cosine edges).
   **Validated against GNPS2's own `ComponentIndex`**: 0 edges span our
   components and the partition agrees 1:1 in both directions.
2. Restricted to features surviving the artifact filter and present in the
   analysis matrix, so families are described in the same feature universe
   as every other tier.
3. Per component: size, per-species blank-clearing member count (paired
   plate rule, union over stages), SIRIUS class composition, example
   structures.
4. Fisher exact (BH-corrected) per component: is its blank-clearing fraction
   above the background rate?
5. Edge `DeltaMZ` classified against common homologous steps (CH2, 2×CH2,
   H2, H2O, C2H4, O, NH, CH3, CO, C2H2) at the workflow's 0.05 Da tolerance.

## Results (2026-09-02)
- Raw network: 9,600 edges / 4,423 nodes / **659 components**, largest 97.
- Artifact-filtered analysis matrix: **551 components**, largest **47**, 226
  with ≥3 members. Families were substantially inflated by isotope/adduct/ISF
  copies of the same molecule cosine-matching itself.
- Background blank-clearing rate: Bd 4.16%, Bsal 7.28%.
- **Blank-enriched components (q<0.05, ≥3 members): 5 (Bd), 8 (Bsal).**
- **Entirely blank-clearing components (≥3 members): 0 (Bd), 2 (Bsal).**

**This is a negative result for the secretion hypothesis** — there is no
population of clean fungal-derived homologous series here.

The species split independently reproduces the peptide finding:

| species | dominant class among blank-enriched components | reading |
|---|---|---|
| Bd | Fatty acids (3/5) — glycerophospholipids, phosphoethanolamine esters | more plausibly cell-derived carryover than secreted |
| Bsal | Amino acids and Peptides (6/8) — `H-Val-Val-Pro-Pro-Phe-OH`, `H-Ala-Pro-Glu-Ala-Val-OH`, `H-Pro-Ser-Pro-Ser-Pro-Ser-al` | proline-rich peptides, matching the composition test via an independent route (network topology, not residue counting) |

DeltaMZ ladders are real but a minority: **237 of 2,781 intra-matrix edges
(8.5%)** classify — CH2 68, 2×CH2 49, H2 49, O 44, C2H2 13, H2O 11, NH 3
(median cosine 0.84–0.89 throughout).

## Caveats
- Components come from the GNPS2 run's own cosine threshold; nothing is
  re-networked here.
- Cosine edges connect features, not guaranteed analogs. The artifact filter
  removes the systematic classes (isotopes, non-default adducts, ISF) but
  co-isolation chimeras remain.
- `all_blank_clearing` is evidence about a family, not proof any member is a
  secreted secondary metabolite: a family of medium-protein digest peptides
  can be entirely blank-clearing if the fungus generates all of them.

## Key outputs
`components.tsv`, `delta_mz_ladders.tsv`, `component_summary.md`,
`component_sizes.{png,pdf}`. Repro: `pixi run network-components`.
