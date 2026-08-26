# Genome-Bioactivity Linkage — Follow-up Summary (2026-08-26)

Autonomous execution of the five follow-up items from the
genome-bioactivity-linkage precision review. All five are complete,
verified, and documented. Nothing committed to git (uncommitted work for
review).

## What was done

### 1. RNA-seq gene counts fixed → expression evidence folded into the linkage

- **Root cause of empty STAR counts:** NCBI `genomic.gff` (GFF3) `exon`
  lines carry `Parent=rna-*` + `locus_tag` but **no `gene_id` attribute**, so
  STAR `--quantMode GeneCounts` assigned every read to `MissingGeneID`.
- **Fix:** `gffread -T` GFF3→GTF conversion (threads gene/mRNA Parent chain
  into `gene_id`/`transcript_id`), then `featureCounts -s 0` (unstranded —
  best-assigning mode by ~2×) over the existing BAMs. **No re-alignment
  needed** (jobs `27777024`, `27778101`).
- Also found & fixed a silent glob bug that had dropped 3 of 8 samples
  (`SRR1301211*`/`SRR2768388*` are not clean SRR prefix families) — BAM paths
  now come from `samples.tsv` (`run_featurecounts.sh`).
- Candidate tables now carry `reference_protein_id`, `ref_locus`,
  `gene_total_raw`, `n_rep_ge_min`, `rna_is_expressed`, `rna_no_evidence`
  (via `build_expression_evidence.py`, `pixi run rna-expression-evidence`).
- Expression summary: Bd 2,634/2,634 candidates' genes detected (Bd baseline
  95% → weak discriminator); Bsal 7,222/8,322 (baseline 63% → informative).
  Tier-1 NRPS genes transcribed both species.

### 2. Tier-1 NRPS genes characterized

- **Bd `FCC698BD_004035-T1`** = BDV3_005439 = NCBI XJO74646.1/XJO74647.1:
  type-I NRPS, C–A(**Ala**, NRPSPredictor)–PCP module + downstream
  transaminase (BDV3_005441). MIBiG best hits: **hassallidin C/D**
  (BGC0000369, NRPS+saccharide), **anachelin** (BGC0002532, NRPS+PKS
  siderophore), gramicidin LgrC-type. ~30% identity → product not
  sequence-predictable.
- **Bsal `F61BA062_001377-T1`** = BSLG_000866 = NCBI KAJ1345353.1: C–A(X)–
  PCP module + flanking large NRPS ORFs (BSLG_000863/867/868) + PPIase.
  MIBiG best hit **nostophycin** (BGC0001029), cyanopeptolin-like simC
  (BGC0000334). Product class not confident.
- Cross-species: both top hits are peptide/siderophore-type NRPS — the same
  ontology family as the dominant "Amino acids and Peptides" class in both
  candidate tables.

### 3. Bsal over-prediction quantified

- DIAMOND all-vs-all self-blastp on both BFD proteomes
  (`results/duplication_check/`, job `27778128`):
  - Bd: 1,234 near-dups (14.7%), 966 strict (11.5%) of 8,396 proteins.
  - Bsal: **5,983 near-dups (30.7%)**, 3,453 strict (17.8%) of 19,449.
  - Bsal also 29.5% of proteins <200 aa (Bd 18.4%).
- BFD emits exactly one transcript per locus (no isoform inflation); Bsal =
  19,449 BFD loci vs 10,867 NCBI genes.
- DeepTMHMM-crasher `F61BA062_016014-T1` (4,777 aa) is **real in both
  annotations** (BSLG_008696 = KAJ1332392.1, 4,681 aa, "hypothetical
  protein"), has no PFAM domains and no internal repeat structure, and is
  transcribed — a genuine uncharacterized giant single-exon ORF, not a BFD
  artifact.

### 4. `parse_pfam_domains.py` offsets verified

- Column offsets correct against real `hmmsearch --domtblout` files
  (`results/pfam_hmmscan/<species>/FCC698BD.domtblout.gz`): fields 0/3/4/6/7
  (protein id, Pfam name, Pfam accession, E-value, score) all land correctly
  for the hmmsearch format the local fallback (and BFD's PFAM module) uses.

### 5. Metabolomics-USI live-spectrum verification scaffold

- Per-species ranked lists of liq-vs-spore-enriched significant features
  joined to SIRIUS identity (`sirius_structure_name/formula/npc_*`,
  `bioactive`, `annotation_origin`) and GNPS library hits
  (`NAME`, `cosine`, `SPECTRUMID`), each with live USI render/mirror/json
  link-outs (`differential_features_primary/liq_enriched_curation/`).
- Validated live: the resolver returns real peak arrays for top hits
  (Bd feature 473 @ m/z 1058.58 → 623 peaks; Bsal 2953 → 481 peaks).
- Stats: Bd 13,170 liq-enriched features (1,659 SIRIUS-structure, 82 lib
  hits); Bsal 7,179 (656 structure, 44 lib hits).
- Top hits are NRPS-peptide-type (e.g. Bd `H-Leu-Leu-Phe-Gly-Nle-Pro-Val-`,
  Bsal `cyclo[Phe-D-Phe-Val-DL-Pro-Pro...]`) — prime MS² verification
  targets vs the tier-1 NRPS BGCs.

## Key artifacts

| Path | Contents |
|------|----------|
| `analysis/rnaseq_expression/RNASEQ_EXPRESSION.md` | gene-count pipeline, bug, strandedness, results, caveats |
| `analysis/rnaseq_expression/scripts/{run_featurecounts.sh, build_expression_evidence.py}` | reproducible counts + expression fold-in |
| `analysis/genome_bioactivity_linkage/results/TIER1_NRPS_CHARACTERIZATION.md` | the two NRPS BGC characterizations (results/ is gitignored; force-add if you want it committed) |
| `analysis/genome_bioactivity_linkage/results/duplication_check/<species>/self.blastp.tsv` | DIAMOND self-blast tables |
| `analysis/genome_bioactivity_linkage/results/{dendrobatidis,salamandrivorans}_candidate_table.tsv` | candidate tables + 4th expression boolean (tracked) |
| `analysis/differential_features_primary/liq_enriched_curation/{dendrobatidis,salamandrivorans}_liq_enriched_top.tsv|_usi.html` | USI curation scaffold + live link grid |
| `analysis/differential_features_primary/scripts/build_usi_curation.py` | USI scaffold generator |
| `pixi.toml` | + `rna-expression-evidence`, `usi-curation` tasks |
| `.living/findings/bsal-overprediction-and-tier1-nrps.md` | F-004 finding |
| `.living/learnings.md`, `.living/log/2026-08-26-001-batrachochytrium-ms.md` | 3 new learnings + session log |

## Repro chain

```bash
pixi run gbl-build-tables && pixi run rna-expression-evidence   # candidate tables + expression
sbatch analysis/rnaseq_expression/scripts/run_featurecounts.sh   # after gffread GTF exists
sbatch analysis/genome_bioactivity_linkage/scripts/run_duplication_check.sh
pixi run usi-curation
```

## Notes / caveats

- `analysis/genome_bioactivity_linkage/results/` and
  `analysis/rnaseq_expression/results/` are **gitignored** (regenerable);
  only the candidate table TSVs are force-tracked. The tier-1
  characterization doc sits under the ignored results tree — force-add if
  it should be committed.
- Two stale tests were updated to match the documented
  `require_extracellular=False` default and the on-disk fallback files;
  `pixi run pytest analysis/genome_bioactivity_linkage/tests` → **49 passed**.
- RNA-seq is presence/absence only (public runs, not the LC-MS growth
  condition); absence is weaker evidence than presence.
- Bsal cross-ref weaker than Bd (raw GenBank vs curated FungiDB).
