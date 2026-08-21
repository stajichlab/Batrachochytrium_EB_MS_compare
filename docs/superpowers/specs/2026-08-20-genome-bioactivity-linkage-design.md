# Genome-to-bioactive-compound linkage — design spec

Status: approved by user, pending Opus bioinformatics review of shell/SLURM
vs. Nextflow implementation choice before writing-plans.

## Purpose

Link SIRIUS/CANOPUS-annotated bioactive metabolite classes (terpenoid,
polyketide, peptide-like, etc.) found in the Bd/Bsal LC-MS/MS data to
candidate biosynthetic genes in the Bd/Bsal genomes — BGCs where they exist
(expected to be sparse/absent in chytrids), dispersed terpene
synthase/PKS/NRPS domains, and secreted-protein candidates for compounds
enriched in the liquid (media) fraction.

## Reference genomes

- **Bd**: JEL423 — BFD `genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/` (8,396 proteins, `.gbk` + `.proteins.fa` already exist from BFD gene prediction; no functional annotation run yet).
- **Bsal**: AMFP13 — BFD `genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/` (19,449 proteins + `.gbk`). Protein count is high for a chytrid (typical ~8–9k genes) — flag as a caveat (possible isoform/fragment inflation) rather than trusting silently; sanity-check against BUSCO completeness stats already in `results/genome_stats_by_name/Batrachochytrium/`.

## Cross-reference source (asymmetric by species)

- **Bd**: JEL423 is the true Bd reference strain (per user; supersedes the older assumption that JAM81 is the FungiDB/RefSeq reference) — so FungiDB/NCBI JEL423 annotation should map directly onto the BFD JEL423 gene models without a separate ortholog-mapping step. Confirm at implementation time which NCBI JEL423 assembly (`GCA_048537975.1_CMM_BatrDend_JEL423_V3` per BFD's `samples.csv`) and which FungiDB release carry curated InterPro/PFAM/GO, and reconcile gene-model coordinates only if the FungiDB/NCBI JEL423 assembly version differs from the one BFD used for gene prediction (fall back to reciprocal-best-hit mapping only if coordinates don't line up).
- **Bsal**: no RefSeq/FungiDB entry exists. Cross-check falls back to the raw NCBI GenBank annotation for `GCA_002006685.2` (AMFP13) instead of a curated source. This is materially weaker than the Bd cross-check — call this out in the writeup, don't present it as equivalent rigor.

## Stage 1 — BGC + dispersed domain detection (per species, on BFD `.gbk`)

- `antismash --taxon fungi --genefinding-tool none --fullhmmer --clusterhmmer --cb-general --pfam2go -c 8 <gbk>` — same invocation as BFD's existing `ANTISMASH_RUN` nextflow module (`~/projects/BFD/Fungi_BFD_runs/nextflow/modules/funannotate/function/ANTISMASH_RUN/main.nf`). Expect few/no clusters called (prior belief, stated by user); a null result here is itself informative and must be reported, not treated as pipeline failure.
- Independent of cluster calls: mine the `--fullhmmer` PFAM domain hits (or a fresh `hmmscan -E 1e-5 Pfam-A.hmm` if antiSMASH's own sweep proves insufficient) for terpene synthase (PF01397, PF03936, PF00494), PKS (PF00109 KS, PF08659, PF00698 AT), NRPS (PF00668 C-domain, PF00501 AMP-binding, PF00550 PCP), DMATS prenyltransferase, and adjacent CYP450 domains — because chytrid secondary metabolism is expected to be dispersed/unclustered rather than organized into classic BGCs.

## Stage 2 — secretion prediction

- **SignalP 6** (HPCC module `signalp/6`) for signal-peptide calls.
- **DeepTMHMM** for transmembrane helix calls, replacing plain TMHMM — port the invocation from `~/projects/nf/nf_funannotate1/modules/local/deeptmhmm_annotation.nf` and `tests/test_deeptmhmm_gpu.sh`: singularity image `/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif`, run as `apptainer exec --nv -B <project_dir> "$SIF" bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta <proteins.fa> --output-dir <outdir>"`, output `TMRs.gff3`. No HPCC Lmod module exists for this — it only runs containerized, and needs a GPU node (`preempt_gpu` partition, not plain `preempt`).
- **PredGPI** for GPI-anchor calls (mirrors the BFD `function/predgpi` category, not yet run for Batrachochytrium).
- Combine: "predicted extracellular" = SignalP+, no TM helix beyond the cleaved signal peptide, no GPI anchor.

## Stage 3 — cross-reference reconciliation

- Bd: since JEL423 is the true reference strain, compare BFD-recompute domain calls directly against FungiDB/NCBI JEL423 curated annotation (coordinate-matched if the assembly version matches; reciprocal-best-hit mapping only as a fallback if it doesn't). Flag domain calls found in both sources as higher-confidence; calls found only in one side are reported but flagged lower-confidence.
- Bsal: same mechanic but against raw AMFP13 GenBank annotation only (no curated FungiDB/InterPro source) — weaker cross-check, noted as such.

## Stage 4 — linking to metabolomics

- Source tables: `analysis/sirius_annotation/sirius_annotations.tsv` (CANOPUS/NPClassifier compound class per feature) joined to `analysis/differential_features_primary/` liquid-vs-spore contrasts, restricted to features enriched in the liquid (media) fraction.
- Compound class → candidate domain family mapping: terpenoid → terpene synthase hits; polyketide → PKS hits; peptide-like/alkaloid-with-amide-bond → NRPS/RiPP-adjacent hits; unmapped classes are left out of the candidate table rather than force-fit.
- Per species/compound-class candidate table, secreted proteins (Stage 2) with matching domain evidence (Stage 1) are ranked by a **composite proxy score** — NOT real co-expression (no RNA-seq/proteomic quantification exists for Bd/Bsal life stages in this project, only a Trinity assembly used as gene-prediction evidence). Score inputs: (a) number/quality of matching biosynthetic domains, (b) SignalP confidence, (c) BGC cluster membership if any (Stage 1), (d) the linked metabolite's differential-abundance significance/fold-change (liq vs spore) as a prioritization weight. Document the score formula plainly as a heuristic ranking, not a causal or expression-based claim.

## Compute

- All SLURM jobs (antiSMASH, SignalP, hmmscan, ortholog BLAST/DIAMOND) submit through `-p preempt -A preempt`.
- DeepTMHMM specifically requires `-p preempt_gpu -A preempt --gres=gpu:1` (GPU + singularity container, no CPU-only path).
- Everything else (parsing, scoring, table-building) runs in this repo's existing pixi env, matching the pattern already used by `analysis/*/scripts/`.

## Output

New `analysis/genome_bioactivity_linkage/` directory:
- `scripts/` — one script per stage, following the existing `analysis/<topic>/scripts/` convention in this repo.
- Per-species domain/secretion/BGC tables (intermediate).
- Final linkage table(s): compound class × candidate gene × evidence/score, one per species.
- `GENOME_BIOACTIVITY_LINKAGE.md` writeup following the `analysis/<topic>/<TOPIC>.md` convention already used by `SIRIUS_ANNOTATION.md`, including the caveats above (Bsal weaker cross-check, proxy score is not real co-expression, Bsal protein-count anomaly).

## Open implementation question (pending review)

Whether Stages 1–4 should be plain shell/SLURM scripts (matching this
repo's existing `analysis/*/scripts/*.py` + sbatch pattern) or whether the
per-genome/per-tool fan-out (antiSMASH × 2 genomes, SignalP × 2, DeepTMHMM ×
2, hmmscan × 2, ortholog BLAST × 1) justifies standing up a small Nextflow
pipeline (this project has no existing Nextflow use, unlike sibling BFD/
nf_funannotate1 projects). A bioinformatics-focused review will assess this
before an implementation plan is written.
