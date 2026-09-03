<!-- BEGIN QUICK REFERENCE -->
# .living/ Index
Last audit: 2026-09-02

| File | Entries | Last updated | Key topics |
|------|---------|--------------|------------|
| conventions.md | 0 sections | 2026-08-19 | — |
| decisions.md | 10 entries | 2026-09-02 | Port EB analysis scripts to the Everything-Bagel feature table with a schema-only adapter, 2026-08-19 — Transfer SIRIUS annotations from the EB project instead of re-running SIRIUS, 2026-08-19 — Keep one row per local feature id in the accumulated annotation table, 2026-08-19 — Native SIRIUS run: charge-1+ targets only, small per-shard jobs, pilot-first, 2026-08-20 — Collapse Sporangium+Mature into a "Developed" stage_group for the primary analysis tier |
| learnings.md | 25 entries | 2026-09-02 | 2026-08-19 — Cross-project SIRIUS annotation transfer for the Everything-Bagel features, 2026-08-19 — Everything-Bagel merges isobaric/co-eluting EB features into single features, 2026-08-19 — Everything-Bagel feature MGF has degenerate blocks (CHARGE=0 / PEPMASS=0.0) even for has_ms2=True features, 2026-08-19 — Native SIRIUS target set reduces mainly on charge state, not sample presence, 2026-08-19 — Everything-Bagel aligned_features.csv area columns differ from the EB quant schema |
| log/ | 8 sessions | 2026-09-02 | batrachochytrium-ms (6), batrachochytrium-eb-ms-compare (2) |
| findings/ | 4 findings across 4 topics | 2026-09-02 | molecular-family-evidence-negative, lifestage-signal-in-spore-fraction, matrix-dominates-bagel-metabolome, bsal-overprediction-and-tier1-nrps |

## Local skills
See `.living/skills/` for project-specific skill packs.
<!-- END QUICK REFERENCE -->

<!-- BEGIN KNOWLEDGE SUMMARY -->
Last summarized: 2026-09-02 (heuristic)

## Tag clusters

- **sirius** (13 entries) — D-2, D-4, D-6, D-7, D-8
- **everything-bagel** (10 entries) — L-9, L-10, L-20, L-21, D-1
- **metabolomics** (10 entries) — D-2, D-3, D-4, D-6, D-7
- **decision** (7 entries) — D-4, D-5, D-6, D-7, D-10
- **genome-bioactivity-linkage** (6 entries) — L-17, L-19, D-8, D-9, D-10
- **analytical-design** (5 entries) — L-7, L-16, L-21, D-5, D-10

## Most recent (10)

- [2026-09-02] L-25: Mycelium's Stop hook and init_repo require Python 3.11 (`datetime.UTC`) but invoke bare `python3`; on a 3.10 box Stop blocks with a message that hides the real error
- [2026-09-02] L-24: The blank-clearing "secreted" peptides are proline-rich at casein/gelatin levels: the secreted metabolome is dominated by medium proteolysis
- [2026-09-02] L-23: Media blanks were 50% of every liq group; and at n=5 "n_significant" is a step function of the feature universe, not an effect size
- [2026-09-02] L-22: Media-blank term applied: 91% of "secreted candidates" were medium; ordinal trend recovers 3-6x more life-stage signal than the collapse
- [2026-09-02] L-21: F-003's "Sporangium-vs-Mature is always 0-significant" is wrong for Bd's spore fraction: 5,507 features separate them
- [2026-09-02] L-20: The media-blank filter exists but is not applied to the secreted-candidate or USI-curation paths; the named MS² priority targets are media peptides
- [2026-08-26] L-19: Bsal BFD proteome over-prediction quantified: 30.7% self-BLAST near-duplicates (Bd 14.7%) + heavy short-protein tail; the DeepTMHMM outlier is real in NCBI too
- [2026-08-26] L-18: SRR id prefixes are not clean glob families: a sloppy ls pattern silently dropped 3 of 8 samples from featureCounts
- [2026-08-26] L-17: STAR --quantMode GeneCounts produced zero gene counts: NCBI GFF3 exons carry no gene_id; fix is gffread GTF + featureCounts on existing BAMs
- [2026-08-25] D-10: Relax genome-bioactivity-linkage's `is_extracellular` gate from a hard filter to an informational column

## By tag

- `sirius`: L-1, L-3, L-4, L-9, L-10, L-13, L-14, L-20, D-2, D-4, D-6, D-7, D-8
- `everything-bagel`: L-3, L-5, L-6, L-7, L-8, L-9, L-10, L-20, L-21, D-1
- `metabolomics`: L-1, L-2, L-3, L-4, L-20, D-2, D-3, D-4, D-6, D-7
- `decision`: D-2, D-3, D-4, D-5, D-6, D-7, D-10
- `genome-bioactivity-linkage`: L-16, L-17, L-19, D-8, D-9, D-10
- `analytical-design`: L-7, L-16, L-21, D-5, D-10
- `differential-abundance`: L-6, L-7, L-21, D-1, D-5
- `collapse`: L-7, L-8, L-21, D-5
- `life-stage`: L-7, L-8, L-21, D-5
- `annotation`: L-9, D-2, D-3
- `hpc`: L-13, L-14, L-15
- `porting`: L-5, L-10, D-1
- `secreted-compounds`: L-7, L-20, D-6
- `slurm`: L-12, L-14, L-15
- `bfd`: D-8, D-9
- `data-engineering`: L-5, L-18
- `data-quality`: L-2, L-3
- `deferral`: L-13, D-7
- `feature-tables`: L-10, L-11
- `featurecounts`: L-17, L-18
- `gpu`: L-15, D-9
- `native-run`: D-4, D-7
- `ordination`: L-6, D-1
- `pipeline-status`: D-8, D-9
- `queue-management`: L-13, D-7
- `secretion-prediction`: L-16, D-10
- `signalp`: L-15, D-9
- `spore-fraction`: L-7, L-21
- `timeout`: L-14, L-15
- `analytical-concordance`: L-6
- `annotation-join`: D-6
- `annotation-transfer`: L-1
- `array-job`: L-14
- `background-subtraction`: L-20
- `bioactivity`: D-6
- `bioinformatics`: L-15
- `blocker`: D-8
- `bsal`: L-19
- `candidate-recovery`: D-10
- `caveat`: L-16
- `chunking`: L-11
- `concurrency`: L-13
- `conflict`: L-2
- `correction`: L-21
- `data-model`: D-3
- `data-schema`: L-8
- `deeptmhmm`: D-8
- `dendrobatidis`: L-21
- `diamond`: L-19
- `duplication`: L-19
- `eb-comparison`: L-6
- `false-positive`: L-20
- `feature-collapse`: L-2
- `feature-matching`: L-1
- `feature-selection`: L-4
- `feature-table`: L-5
- `file-size-limit`: L-11
- `gene-content`: L-19
- `gene-counts`: L-17
- `gff3`: L-17
- `gffread`: L-17
- `git-push`: L-11
- `github`: L-11
- `glob`: L-18
- `gnps`: L-1
- `gtf`: L-17
- `hpcc-short`: L-13
- `html`: L-10
- `infrastructure`: L-9
- `interactive`: L-10
- `is_extracellular`: D-10
- `isobaric`: L-2
- `login-token`: L-13
- `media-blank`: L-20
- `mgf`: L-3
- `monitoring`: L-12
- `ms2-cosine`: L-1
- `mycelium`: L-12
- `native-merge`: L-9
- `over-prediction`: L-19
- `pfam`: L-16
- `pilot`: L-9
- `power`: D-5
- `reproducibility`: D-1
- `retry`: L-14
- `rnaseq`: L-17
- `rollup`: L-11
- `sample-metadata`: L-8
- `sample-sheet`: L-18
- `schema`: L-5
- `scoping`: L-4
- `session-lifecycle`: L-12
- `sharding`: D-4
- `silent-failure`: L-18
- `srr`: L-18
- `standalone-repo`: D-1
- `star`: L-17
- `stop-hook`: L-12
- `tooling-gotcha`: L-12
- `unblocked`: D-9
- `usi-curation`: L-20
- `validation`: L-3

_Heuristic clustering: tags with ≥2 entries, top 6 by count. To fetch matching entries: `python3 "$(cat .mycelium/plugin-root)/skills/core/scripts/recall_lessons.py" --living-dir <path> --tag <tag>` or `--id L-N`._
<!-- END KNOWLEDGE SUMMARY -->
