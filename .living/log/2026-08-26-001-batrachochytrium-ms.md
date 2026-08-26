---
session_id: 2026-08-26-001
project: batrachochytrium-ms
branch: "main"
started: 2026-08-26T00:10:00-0700
ended:
duration_minutes:
files_changed:
---

## Session Log

### Resume from 2026-08-25-003; autononomous follow-up on genome-bioactivity linkage precision

User asked to execute the 5 follow-up items from the genome-bioactivity
linkage review autonomously.

1. **RNA-seq gene counts fixed** — root cause: STAR `--quantMode GeneCounts`
   dropped all gene counts because NCBI GFF3 exon lines have no `gene_id`
   attribute. Fixed with `gffread -T` GTF conversion + `featureCounts -s 0`
   on existing BAMs (job 27777024). First run's glob silently dropped 3/8
   samples (SRR prefix family bug) — found & fixed (job 27778101). Folded
   `rna_is_expressed` / `gene_total_raw` / `n_rep_ge_min` / `rna_no_evidence`
   into both candidate tables via `build_expression_evidence.py`. Bd
   2,634/2,634 expressed (weak — Bd baseline 95%); Bsal 7,222/8,322
   (informative — baseline 63%). Tier-1 NRPS genes transcribed both species.
2. **Tier-1 NRPS characterization** — Bd FCC698BD_004035-T1 = BDV3_005439 =
   XJO74646.1/47.1, C-A(**Ala**)-PCP NRPS + downstream transaminase, MIBiG
   hassallidin C/D (BGC0000369) / anachelin (BGC0002532); Bsal
   F61BA062_001377-T1 = BSLG_000866 = KAJ1345353.1, C-A(X)-PCP + flanking
   NRPS ORFs + PPIase, MIBiG nostophycin (BGC0001029) / cyanopeptolin-like
   (BGC0000334). Written to `results/TIER1_NRPS_CHARACTERIZATION.md`.
3. **Bsal duplication check** (DIAMOND self-blast, job 27778128) — Bsal
   30.7% near-dups vs Bd 14.7%; strict 17.8% vs 11.5%; single transcript per
   BFD locus; outlier F61BA062_016014-T1 real in NCBI too (BSLG_008696,
   4,681 aa), no PFAM/repeats, transcribed.
4. **parse_pfam_domains.py offsets verified** against real hmmsearch
   domtblout — correct (fields 0/3/4/6/7).
5. **metabolomics-USI curation scaffold** — per-species ranked liq-enriched
   top features with live USI render/mirror/json links; endpoints validated
   (real peak arrays). Bd 13,170 features (82 lib hits), Bsal 7,179 (44).

Docs updated: RNASEQ_EXPRESSION.md (new), TIER1_NRPS_CHARACTERIZATION.md
(new), GENOME_BIOACTIVITY_LINKAGE.md (Follow-ups resolved), ANALYSIS_MANIFEST.md
(3 entries), pixi.toml (2 tasks), F-004 finding + FINDINGS_REGISTRY,
learnings (3 entries).

### Files Modified (selection)
- analysis/rnaseq_expression/scripts/{run_featurecounts.sh,build_expression_evidence.py}
- analysis/genome_bioactivity_linkage/results/{dendrobatidis,salamandrivorans}_candidate_table.tsv
- analysis/genome_bioactivity_linkage/results/TIER1_NRPS_CHARACTERIZATION.md
- analysis/genome_bioactivity_linkage/results/duplication_check/*
- analysis/differential_features_primary/liq_enriched_curation/*
- analysis/rnaseq_expression/RNASEQ_EXPRESSION.md
- analysis/GENOME_BIOACTIVITY_LINKAGE.md, ANALYSIS_MANIFEST.md, pixi.toml
- .living/findings/bsal-overprediction-and-tier1-nrps.md, FINDINGS_REGISTRY.md, learnings.md
