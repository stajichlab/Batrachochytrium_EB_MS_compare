# Genome-Bioactivity Linkage

## Purpose

Link mass-spec metabolite classes (from SIRIUS/CANOPUS compound-class
annotation, see `analysis/sirius_annotation/SIRIUS_ANNOTATION.md`) to
candidate biosynthetic genes in the two chytrid genomes analyzed here —
*Batrachochytrium dendrobatidis* JEL423 and *B. salamandrivorans* AMFP13 —
producing a tiered, per-species candidate-gene table for each
background-filtered, class-mapped liquid-fraction compound.

The pipeline runs on top of a separate **BFD (Fungal Genome Annotation)**
run (`/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs`) that produces
gene predictions, PFAM domain scans, antiSMASH BGC calls, and secretion
predictors (SignalP/DeepTMHMM/PredGPI) for both genomes. This project's
`analysis/genome_bioactivity_linkage/` code is glue: it parses BFD's Stage
1–2 outputs, cross-references them against a curated NCBI reference
annotation (Stage 3), and joins the result against this project's
background-filtered, SIRIUS-annotated, differentially-abundant liquid
compounds (Stage 4) to rank candidate genes per compound.

## Status

**Implemented and unit-tested; not yet run against real BFD output.**

All ten build tasks in the implementation plan are complete:

- Task 1: shared path/domain constants (`paths.py`, `domain_families.py`)
- Task 2: media-companion background filter (`background_subtraction.py`)
- Task 3: DeepTMHMM wrapper + parser (`run_deeptmhmm.sh`, `parse_deeptmhmm.py`)
- Task 4: NCBI reference-annotation fetch (`fetch_reference_annotation.sh`) —
  **already run**, real data on disk under `results/reference_annotation/`
- Task 5–6: DIAMOND RBH ortholog mapping (`run_rbh.sh`, `parse_rbh.py`)
- Task 6–7: BFD PFAM domtblout parser + antiSMASH BGC/fullhmmer JSON parser
  (`parse_pfam_domains.py`, `parse_antismash_clusters.py`)
- Task 8: SignalP/DeepTMHMM/PredGPI secretion merge (`merge_secretion.py`)
- Task 9: tiered compound-to-gene linking (`link_compounds_to_genes.py`)
- Task 10 (this task): end-to-end driver (`build_linkage_tables.py`), pixi
  task registration, this writeup

**Blocked on:** the BFD `--taxon GENUS:Batrachochytrium` functional-annotation
run (PFAM `hmmscan`/`hmmsearch` domtblout, antiSMASH `--fullhmmer` JSON,
SignalP, PredGPI) has not finished producing output for either species as of
this writing. `paths.find_bfd_output` / `paths.bfd_antismash_json` raise
`FileNotFoundError` by design when their target files are absent — this is
correct "skip cleanly / rerun on demand" behavior, not a bug. Once that BFD
run completes:

```bash
pixi run gbl-fetch-reference   # already run; re-runs are a no-op (skip-if-exists)
sbatch analysis/genome_bioactivity_linkage/scripts/run_deeptmhmm.sh
sbatch analysis/genome_bioactivity_linkage/scripts/run_rbh.sh
# after both SLURM jobs complete:
pixi run gbl-build-tables
```

At that point, spot-check that `parse_pfam_domains.py`'s hardcoded domtblout
column offsets and `merge_secretion.py`'s SignalP/PredGPI loader (plain
`pd.read_csv(..., sep="\t")` with the `protein_id, is_signal_peptide,
cleavage_site` / `protein_id, has_gpi_anchor` column contract Task 8
assumed) actually match BFD's real output format, and adjust
`build_linkage_tables.py`'s loader lines if not.

## Datasets

- BFD gene predictions / PFAM domtblout / antiSMASH JSON / SignalP / PredGPI
  — `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/genome_annotation/<out>/`
  and `.../results/function/<kind>/` (see `paths.py`).
- NCBI RefSeq/GenBank reference annotation — fetched by Task 4, on disk at
  `analysis/genome_bioactivity_linkage/results/reference_annotation/<species>/`
  (`GCA_048537975.1_CMM_BatrDend_JEL423_V3` for Bd,
  `GCA_002006685.2_Batr_sala_V2` for Bsal).
- `analysis/sirius_annotation/sirius_annotations.tsv` — SIRIUS/CANOPUS
  compound-class annotations (see that project's `SIRIUS_ANNOTATION.md` for
  its own status/snapshot). **SIRIUS coverage snapshot used for this run:**
  1,885 annotated local features as of 2026-08-20 (1,773 transferred + 112
  native pilot, job `27605104`) — the full native SIRIUS run remains
  deferred; if `sirius_annotations.tsv` is regenerated with a larger native
  run, this note must be re-stated with the new feature count/date, since
  `build_linkage_tables.py` reads it fresh on every invocation.
- `analysis/differential_features_primary/all_significant_features_summary.tsv`
  — 103,638-row primary differential-features table; this pipeline restricts
  to rows whose `comparison` string contains `"liq"` (the liquid-fraction
  contrasts), matching Stage 4's liquid-fraction focus, and keeps the
  lowest-`q_value` row per compound `row_id` as that compound's tie-breaker
  `log2fc`/`q_value`.
- `data/metdata/curated_gnps_metadata.tsv` + this project's aligned feature
  table (`fungal_over_blank_ratio` background filter, Task 2).

## Method

Per species (`dendrobatidis`, `salamandrivorans`):

1. **PFAM domain classification** (Stage 1): parse BFD's `hmmscan`/
   `hmmsearch --cut_ga` domtblout via `parse_pfam_domains.parse_domtblout`,
   classify each hit into a biosynthetic domain family
   (`terpene_synthase`, `squalene_phytoene_synthase`, `pks`, `nrps`,
   `dmats_prenyltransferase`, `p450`) via `domain_families.classify_pfam`.
2. **BGC-region membership** (Stage 1): load antiSMASH's BGC regions
   (`parse_antismash_clusters.load_regions`) and, **per protein**, its real
   genomic coordinates parsed once from the BFD `.gbk` file
   (`SeqIO.parse(..., "genbank")`, `feature.qualifiers["locus_tag"][0]` +
   `feature.location.start`/`.end`, 0-based, matching antiSMASH's own region
   coordinate convention). `has_bgc_context` is `protein_in_bgc(locus_tag,
   (start, end), regions, record_id)` evaluated against those real
   coordinates — **not** membership in antiSMASH's `full_hmmer` hit list.
   antiSMASH's `full_hmmer` module is a genome-wide internal PFAM scan (most
   antiSMASH JSON records have an *empty* `areas` list yet still contribute
   `full_hmmer` hits across the whole genome), so "protein appears in
   `full_hmmer` hits" is not a BGC-region signal and would mislabel nearly
   every domain hit as (false) BGC context if used directly.
3. **Secretion prediction** (Stage 2): merge SignalP + DeepTMHMM (TM-helix-
   outside-signal rule) + PredGPI via `merge_secretion.predicted_extracellular`.
4. **Ortholog cross-reference** (Stage 3): reciprocal-best-hit DIAMOND blastp
   against the NCBI reference proteome (`parse_rbh.reciprocal_best_hits`);
   `is_cross_ref_confirmed` is set for proteins with a confirmed RBH.
5. **Compound-side filtering** (Stage 4): background-filter liquid-fraction
   features (fungal-over-blank ratio ≥ 2×, per life stage) via
   `background_subtraction.fungal_over_blank_ratio`; map each surviving
   compound's SIRIUS NPC pathway/class to a domain family via
   `domain_families.COMPOUND_CLASS_TO_FAMILY`; join each compound `row_id` to
   its most-significant liq-fraction differential-abundance result
   (lowest `q_value` among `comparison` strings containing `"liq"`) from
   `all_significant_features_summary.tsv` — compounds with **no** matching
   liq-fraction contrast row are excluded from the candidate table entirely
   (never shown to be differentially abundant in the liquid fraction, so no
   fabricated tie-breaker value is assigned).
6. **Tiered linking** (Stage 4): `link_compounds_to_genes.build_candidate_table`
   matches each compound's domain family against extracellular
   (`is_extracellular`) proteins of that family, assigns a tier
   (1 = BGC context **and** cross-ref confirmed; 2 = either one; 3 =
   neither), and sorts within a compound by tier, then by
   `abs(compound_log2fc)` descending (real differential-abundance magnitude,
   not a placeholder) — a **tiered/lexicographic ranking**, deliberately
   never collapsed into a single weighted composite score.

Output: `analysis/genome_bioactivity_linkage/results/<species>_candidate_table.tsv`,
one row per (compound, candidate protein) pair.

## Known caveats

1. **Bsal's cross-check is weaker than Bd's.** Bd JEL423 has a curated
   FungiDB/RefSeq annotation to cross-reference against; Bsal AMFP13's NCBI
   entry (`GCA_002006685.2`) is raw GenBank annotation only, with no curated
   FungiDB source. `is_cross_ref_confirmed` for Bsal candidates should be
   read as weaker evidence than the same flag for Bd candidates.
2. **The tiered ranking is a heuristic, not real co-expression evidence.**
   No RNA-seq or proteomic quantification exists for these life stages in
   this project — `has_bgc_context` and `is_cross_ref_confirmed` are
   sequence/genomic-context proxies for "this gene plausibly makes this
   compound class," not a measured expression or co-occurrence signal.
   Tier 1 candidates are the most defensible leads to prioritize for
   follow-up, not confirmed producers.
3. **Bsal protein-count anomaly (unresolved, two independent numbers now
   known).** BFD's own gene-prediction proteome for Bsal AMFP13 has
   **19,449 proteins**, vs Bd JEL423's **8,396** — roughly 2× expected.
   Separately, NCBI's own RefSeq/GenBank annotation for the *same* Bsal
   AMFP13 assembly (`GCA_002006685.2`, fetched by Task 4) independently shows
   **10,867 proteins** vs Bd JEL423's **8,588** — still elevated, but far
   less extreme than the BFD figure. This suggests some genuine biology
   (Bsal may simply have more genes than Bd) plus some degree of
   BFD-specific over-prediction or gene-model duplication. Recommended next
   diagnostic step: a duplication-rate check (e.g. BUSCO duplication rate, or
   an all-vs-all self-blast on the BFD proteome) once the BFD run completes —
   the natural place to run it is Task 5/6's cross-reference step, per the
   original plan's Stage 3.
4. **SIRIUS coverage is a snapshot, not final.** As of this writing,
   `sirius_annotations.tsv` reflects the transfer path plus a 112-feature
   native SIRIUS pilot (1,885 annotated features total; see
   `analysis/sirius_annotation/SIRIUS_ANNOTATION.md`). The full native
   SIRIUS run is deferred. Because `build_linkage_tables.py` reads
   `sirius_annotations.tsv` fresh on every invocation, candidate-table
   contents will change (likely grow) once a larger native run is folded in
   — re-state the snapshot date/feature count here whenever Stage 4 is
   rerun.
5. **Media is a species confounder**, inherited from the underlying
   metabolomics data (Bd `1% Tryptone` vs Bsal `50% TGHL`) — see the
   project's top-level `CLAUDE.md`. The background filter compares each
   species against its own matched `C_liq` companion blank, so this does not
   bias `passes_background_filter`, but cross-species compound-class
   comparisons should still be read with the media difference in mind.
