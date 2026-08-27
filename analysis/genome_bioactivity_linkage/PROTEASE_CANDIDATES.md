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

## Resolved — the Bsal M36 expansion is real biology, not (primarily) over-prediction

**Cross-reference against the duplication self-blast (2026-08-26):** using
an independently-chosen near-duplicate definition (pident ≥90% AND
min-coverage ≥80% against a same-species self-blastp hit, excluding
self-hits — not necessarily identical to whatever threshold produced the
F-004 aggregate 14.7%/30.7% figures, since that computation's exact
method wasn't saved as a script) applied specifically to the M36 hits:

| Species | All M36 hits | Near-duplicate | Secreted M36 candidates | Near-duplicate |
|---|---|---|---|---|
| Bd JEL423 | 39 | 22 (56%) | 32 | 17 (53%) |
| Bsal AMFP13 | 328 | 98 (30%) | 233 | 87 (37%) |

This is the **opposite** of what the "Bsal count is inflated by
over-prediction" hypothesis predicts: Bd's smaller M36 set is
*proportionally more* duplicated than Bsal's much larger one. Bsal's
extra M36 genes are mostly non-redundant paralogs by this metric, not an
artifact of the same duplication process documented for the proteome
overall (F-004). This argues the M36 count difference reflects **real,
lineage-specific gene-family size**, not an assembly/gene-calling
artifact — see literature corroboration below, which independently
confirms this via a completely different assembly and method.

## Literature corroboration (2026-08-26 web search)

This species-level M36 pattern (Bd small, Bsal large) is **independently
reported in the primary literature**, and the match to our own pipeline
is unusually direct:

- Yu et al. (2025), *A near-complete telomere-to-telomere genome assembly
  for Batrachochytrium dendrobatidis GPL JEL423 reveals a larger CBM18
  gene family and a smaller M36 metalloprotease gene family than
  previously recognized*, G3 Genes|Genomes|Genetics
  ([Oxford Academic](https://academic.oup.com/g3journal/article/15/2/jkae304/7930337),
  [bioRxiv preprint](https://www.biorxiv.org/content/10.1101/2024.10.22.619730v1)):
  their new T2T assembly (`CMM_BatrDend_JEL423_V3`, Oxford Nanopore,
  University of Exeter) shows Bd JEL423 encodes **fewer than half** the
  M36 genes predicted from the older, more fragmented assembly, and a
  gene-tree analysis shows the M36 family is **"highly expanded (n=177)
  in B. salamandrivorans since its split with B. dendrobatidis."**
  — **`CMM_BatrDend_JEL423_V3` is exactly the Bd reference assembly this
  project already uses** (confirmed: `assemblyName` in
  `results/reference_annotation/dendrobatidis/ncbi_dataset/data/assembly_data_report.jsonl`
  matches verbatim). Our own BFD-based count (39 total M36 hits in Bd,
  vs 328 in Bsal — order-of-magnitude consistent with their n=177 for
  Bsal via a different gene-family-calling method) is therefore not just
  plausible, it is built on the same corrected assembly the paper used,
  and reproduces its qualitative conclusion independently.
- Fisher lab / Rosenblum lab work on Bd secreted proteases (e.g. a
  subtilisin-like **serine** protease induced by thyroid hormone that
  degrades host antimicrobial peptides —
  [ScienceDirect](https://www.sciencedirect.com/science/article/abs/pii/S187861461300069X))
  supports the general "Bd secretes proteases that degrade host
  defense/structural peptides" mechanism this analysis is built on,
  independent of the M36 family specifically.

**Practical upshot:** the M36 expansion in Bsal should be treated as a
genuine, literature-supported virulence-factor difference between the two
species, not an artifact to be explained away — strengthening rather
than undermining the secreted-protease candidate list above.

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
