# Secreted Protease Candidates (MEROPS)

## Purpose

Companion analysis to the main compound-linkage candidate tables
(`GENOME_BIOACTIVITY_LINKAGE.md`) and the Tier-1 NRPS characterization
(`TIER1_NRPS_CHARACTERIZATION.md`). Proteases don't biosynthesize a
small-molecule compound the way NRPS/PKS/terpene-synthase genes do, so
they are **not** joined to SIRIUS compound classes — there is no
"protease produces this metabolite" relationship to encode here. Instead
this ranks BFD-predicted proteins with a MEROPS peptidase-family hit by
the same secretion / orthology / RNA-expression evidence used elsewhere
in this pipeline, on two premises:

1. Fungal pathogens (including other chytrids and dermatophytes)
   routinely use **secreted proteases** as virulence factors, degrading
   host structural/skin proteins.
2. Secreted proteolysis of host/media proteins is a **plausible source**
   of some of the "Amino acids and Peptides"-class liquid-enriched
   metabolomics signal (proteolytic fragments), independent of — and
   competing with — the Tier-1 NRPS biosynthesis hypothesis.

## Method

Per species (`dendrobatidis`, `salamandrivorans`):

1. **MEROPS classification**: `run_merops.sh` runs `blastp` (BFD's own
   invocation, mirrored exactly: `merops_scan.lib`, `-evalue 1e-10
   -max_target_seqs 10 -seg yes -soft_masking true -use_sw_tback`) against
   the BFD-predicted proteome, since BFD's own shared `merops` run has
   zero output for either Batrachochytrium locustag (confirmed
   2026-08-26 — same situation PFAM/SignalP/PredGPI/antiSMASH were in).
   `parse_merops.best_merops_hit` keeps the lowest-e-value hit per
   protein and annotates it with MEROPS clan/family/catalytic-type
   (serine/cysteine/aspartic/metallo/threonine/glutamic/asparagine, or
   `inhibitor` for MEROPS "I"-prefix families, which are peptidase
   *inhibitors*, not peptidases — flagged via `is_inhibitor_family`
   rather than silently dropped).
2. **Secretion prediction**: identical `merge_secretion.predicted_extracellular`
   (SignalP + DeepTMHMM TM-helix-outside-signal rule + PredGPI) used in
   the main pipeline.
3. **Ortholog cross-reference**: identical RBH-to-NCBI-reference
   confirmation used in the main pipeline.
4. **RNA-seq expression**: identical RBH → NCBI locus → `featureCounts`
   join used in `build_expression_evidence.py` (presence/absence only,
   not condition-matched — see that script's caveat).
5. No hard filter is applied — `build_protease_candidates.py` emits every
   MEROPS-hit protein with `is_extracellular`/`is_inhibitor_family` as
   informational columns, sorted `(is_extracellular desc, rna_is_expressed
   desc, evalue asc)`, so nothing is silently excluded.

Output: `results/{dendrobatidis,salamandrivorans}_protease_candidates.tsv`.

## Results (2026-08-26)

| Species | MEROPS-hit proteins | Secreted, non-inhibitor | ...also expressed |
|---|---|---|---|
| Bd JEL423 | 247 | 45 | 20 |
| Bsal AMFP13 | 653 | 247 | 22 |

**The secreted-candidate set is dominated by MEROPS family M36** (the
fungalysin/deuterolysin metalloprotease family — well documented as a
secreted virulence-associated protease family in pathogenic fungi,
including dermatophyte keratinases):

| Species | M36 among secreted candidates | M36 also expressed |
|---|---|---|
| Bd JEL423 | 32 / 45 (71%) | 9 |
| Bsal AMFP13 | 233 / 247 (94%) | 13 |

All top-ranked secreted candidates (lowest e-value) in both species are
RBH-cross-ref-confirmed against the NCBI reference proteome.

## Open question — is the Bsal M36 count real expansion or over-prediction?

Bsal's raw MEROPS-hit and secreted-candidate counts are ~3-7x Bd's, in
the same direction as (and possibly compounding) the already-documented
Bsal protein-count/duplication anomaly (F-004,
`.living/findings/bsal-overprediction-and-tier1-nrps.md`: 30.7%
near-duplicate proteins vs Bd's 14.7%). This has **not** been
disambiguated here — it is plausible that Bsal genuinely has an expanded
fungalysin repertoire (a real, biologically interesting virulence-gene
expansion, consistent with Bsal's more severe/rapid pathogenicity in
some hosts), that near-duplicate gene models are inflating the count, or
both. Cross-referencing the M36-hit protein ids against
`results/duplication_check/salamandrivorans/self.blastp.tsv`'s
near-duplicate pairs would settle this and is the natural next step.

## Caveats

- Same RNA-seq caveat as the main pipeline: expression evidence is
  presence/absence from public runs not grown in the liquid-culture
  condition sampled by the mass-spec data; absence is weaker evidence
  than presence.
- Same Bsal cross-reference caveat as the main pipeline: Bsal's NCBI
  annotation is raw GenBank, not a curated FungiDB source, so
  `is_cross_ref_confirmed` is weaker evidence for Bsal than for Bd.
- This table has not been linked to specific liq-enriched peptide/amino-
  acid metabolomics features the way the Tier-1 NRPS candidates were
  (`liq_enriched_curation/`) — that closure step (does the peptide
  fragment profile look more like proteolytic degradation products than
  an NRPS product?) is a natural follow-up, not yet done.
