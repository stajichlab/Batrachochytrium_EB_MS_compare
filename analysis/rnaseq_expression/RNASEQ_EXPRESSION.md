# RNA-seq Expression Evidence

## Purpose

Add presence/absence expression evidence for the BFD gene models in the
genome-bioactivity-linkage candidate-gene tiering. The original
`genome_bioactivity_linkage` tiering rests on three sequence/genomic-context
booleans (`has_bgc_context`, `is_cross_ref_confirmed`, `is_extracellular`);
RNA-seq gives a **fourth, measured** boolean (`rna_is_expressed`) derived from
public expression data, the practical "more precision" step for tying
candidate biosynthetic genes to observed compounds.

**Caveat carried from design (commit 13f1c83):** the 8 public RNA-seq runs are
NOT from the liquid-culture growth condition sampled by the mass-spec data, so
this is *expression in some lab condition* — usable as presence/absence
("is the locus a real, transcribed gene?") but NOT condition-matched
co-expression with the metabolites. Absence is weaker evidence than presence
because of that mismatch.

## Status

**counts complete, folded into the candidate tables (2026-08-25/26).**

- 8 samples downloaded / aligned: Bd JEL423 3× (SRR27683879-81), Bsal AMFP13
  5× (SRP291769 SRR13012113/117/121/125/129).
- STAR alignments (jobs `27752021` + `27754226` for the 3 Bd tasks that needed
  `--limitBAMsortRAM` against the default estimate) all produced
  coordinate-sorted BAMs.
- **Gnarly bug found & fixed:** STAR `--quantMode GeneCounts` produced **no
  gene counts at all** — every `ReadsPerGene.out.tab` had only the 4 summary
  rows + `MissingGeneID`. Root cause: the NCBI `genomic.gff` (GFF3) `exon`
  lines carry `Parent=rna-*` + `locus_tag` but **no `gene_id` attribute**, and
  STAR's GeneCounts tallies on `gene_id`. Fixed by converting each GFF3 to a
  proper GTF with `gffread -T` (threads the `gene`/`mRNA` Parent chain into
  `gene_id`/`transcript_id`) and recounting the existing BAMs with
  `featureCounts` (`run_featurecounts.sh`, job `27777024`) — no realignment
  needed.
- **Strandedness:** unstranded (`featureCounts -s 0`) assigns ~2× more reads
  than `-s 1/-s 2` for both species, so `-s 0` counts are used.
- Expression folded via `build_expression_evidence.py` (`pixi run
  rna-expression-evidence`): each BFD candidate protein is mapped through its
  reciprocal-best-hit NCBI protein to that protein's NCBI `locus_tag`
  (`protein_id`→`locus_tag` from `genomic.gff` CDS lines), then joined to the
  per-gene fragment counts.

## Datasets / Inputs

- `data/raw/rnaseq_srp291769/<SRR>/{1,2}.fastq.gz` — downloaded SRA runs
  (`download_sra.sh`).
- `analysis/rnaseq_expression/results/star_index/<species>/` — STAR indices
  built from the NCBI reference genome + annotation
  (`build_star_index.sh`); the STAR index GTF source was the NCBI `genomic.gff`
  (see the gene_id bug above for why STAR counts were dropped).
- `analysis/rnaseq_expression/results/star_align/<SRR>/Aligned.sortedByCoord.out.bam`
  — alignments (`align_star.sh`).
- `analysis/rnaseq_expression/results/gtf/<species>.gtf` — gffread-converted
  GTFs (this is what makes counting work).
- NCBI `genomic.gff` (per-species) — protein_id→locus_tag map.
- `analysis/genome_bioactivity_linkage/results/rbh/<species>/{fwd,rev}.tsv` —
  BFD↔NCBI reciprocal-best-hit pairs for the protein-side join.

## Method

1. `gffread -T -o <species>.gtf <species>.genomic.gff` (module `gffread/0.12.7`).
2. `sbatch analysis/rnaseq_expression/scripts/run_featurecounts.sh` (array
   0-1): for each species, `featureCounts -s {0,1,2} -p --countReadPairs -T 8`
   over that species' 3/5 BAMs into
   `results/gene_counts/<species>/counts_s{0,1,2}.txt` + `.summary`.
3. Pick `-s 0` (highest Assigned; this dataset is unstranded).
4. `pixi run rna-expression-evidence` → writes
   `results/gene_counts/<species>/gene_expression.tsv` (per-gene raw counts
   and per-replicate presence per species) and appends to each candidate
   table (`results/<species>_candidate_table.tsv`):
   - `reference_protein_id` — the BFD protein's RBH (else top-hit) NCBI protein
   - `ref_locus` — that protein's NCBI `locus_tag`
   - `gene_total_raw` — summed fragment counts across that species' runs
   - `n_rep_ge_min` — how many replicates have ≥ MIN_COUNTS (10) reads
   - `rna_is_expressed` — `n_rep_ge_min >= MIN_REPS` (default 1)
   - `rna_no_evidence` — candidate whose protein had no NCBI ortholog/locus

## Results

| Species | BFD genes | genes detected ≥10 reads in ≥1 rep | candidate rows | candidates `rna_is_expressed` |
|---------|-----------|------------------------------------|----------------|------------------------------|
| Bd JEL423 | 7,308 gene_ids | 6,954 (95.2%) | 2,634 | 2,634 (100%) |
| Bsal AMFP13 | 10,867 gene_ids | 6,837 (62.9%) | 8,322 | 7,222 (86.8%; 275 rows have no NCBI locus) |

Tier-1 NRPS genes:

- **Bd `FCC698BD_004035-T1`** (BDV3_005439; XJO74646.1): 2,534 raw counts,
  3/3 reps ≥10 → **expressed**.
- **Bsal `F61BA062_001377-T1`** (BSLG_000866; KAJ1345353.1): 739 raw counts,
  5/5 reps ≥10 → **expressed**.
- The DeepTMHMM-excluded outlier `F61BA062_016014-T1` (BSLG_008696;
  KAJ1332392.1): 274/179 counts in the Bsal reps → transcribed, consistent
  with it being a real (if odd) annotated gene, not a truncated artifact.

Interpretation guardrail: Bd's 95% global detection rate makes
`rna_is_expressed` a weak discriminator for Bd (nearly everything is detected
at 20M+ read depth); the column is most informative for Bsal (63% baseline),
where 87% of candidates are measurably transcribed. Treat "expressed" here as
"this locus is a transcribed gene", never as measured upregulation in the
experimental (liquid-culture) condition.

## Known caveats

1. Gene-level detection is **not** condition-matched to the LC-MS growth
   condition; absence is weaker evidence than presence.
2. The mapping BFD protein → NCBI protein relies on RBH/top-hit identity; for
   Bd JEL423 the NCBI annotation is curated (good), for Bsal AMFP13 it is
   raw GenBank (weaker) — same cross-ref weakness as the linkage pipeline.
3. `counts_s0.txt` counts fragments (pairs); longer genes and higher-depth
   libraries (Bd) systemically detect more genes, hence the species
   difference in baseline detection.
