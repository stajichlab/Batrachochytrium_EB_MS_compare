# Genome-to-bioactive-compound linkage — design spec

Status: approved by user; Opus bioinformatics review completed
(2026-08-20) and folded in below. Ready for writing-plans.

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

## Stage 1 — BGC + dispersed domain detection (run inside the BFD pipeline, not this repo)

BFD already has both genomes as rows in `samples.csv`
(`GCA_048537975.1_CMM_BatrDend_JEL423_V3`, `GCA_002006685.2_Batr_sala_V2`)
with `predict_results/` populated and no functional stage run yet, plus
ready-made `ANTISMASH_RUN`, `RUN_PFAM` nextflow modules and a
`--taxon GENUS:Batrachochytrium` filter to scope a run to exactly these two
genomes. Run it there — writing functional annotation into BFD's tree
benefits every downstream BFD user instead of forking a private copy — then
read the outputs back from this repo:

- `antismash --taxon fungi --genefinding-tool none --fullhmmer --clusterhmmer --cb-general --pfam2go -c 8 <gbk>` via BFD's existing `ANTISMASH_RUN` module (`~/projects/BFD/Fungi_BFD_runs/nextflow/modules/funannotate/function/ANTISMASH_RUN/main.nf`), launched with `sbatch nextflow/run_functional.sh --taxon GENUS:Batrachochytrium` (adjust flags as BFD's actual driver script requires at implementation time). Expect few/no clusters called (prior belief, stated by user); a null result here is itself informative and must be reported, not treated as pipeline failure.
- Independent of cluster calls: mine PFAM domain hits via BFD's `RUN_PFAM` module using **`hmmsearch --cut_ga Pfam-A.hmm`** (PFAM's own gathering thresholds — not an arbitrary `-E 1e-5` cutoff) for terpene synthase (PF01397, PF03936), PKS (PF00109 KS, PF08659, PF00698 AT), NRPS (PF00668 C-domain, PF00501 AMP-binding, PF00550 PCP), DMATS prenyltransferase, and adjacent CYP450 domains. **PF00494 is squalene/phytoene synthase, not a terpene cyclase** — track it as its own category, not folded into "terpene synthase"; chytrid TPS genes may simply be absent from PF01397/PF03936 (a true negative, not a pipeline gap).

## Stage 2 — secretion prediction (SignalP/PredGPI via BFD; DeepTMHMM standalone)

- **SignalP 6** and **PredGPI** via BFD's existing `SIGNALP_RUN` / `RUN_PREDGPI` modules, same `--taxon GENUS:Batrachochytrium`-scoped run as Stage 1 (mirrors the BFD `function/signalp`, `function/predgpi` categories, not yet populated for Batrachochytrium).
- **DeepTMHMM** has no BFD/Lmod equivalent yet, so run it standalone in this repo: port the invocation from `~/projects/nf/nf_funannotate1/modules/local/deeptmhmm_annotation.nf` and `tests/test_deeptmhmm_gpu.sh` — singularity image `/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif`, run as `apptainer exec --nv -B <project_dir> "$SIF" bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta <proteins.fa> --output-dir <outdir>"`, output `TMRs.gff3`. Needs a GPU node (`preempt_gpu`, not plain `preempt`).
- Combine: "predicted extracellular" = SignalP+, no GPI anchor, **and no TM helix outside the cleaved signal-peptide region** — define this overlap rule explicitly in the merge script (a TM helix call that falls entirely within the SignalP-cleaved N-terminal region must NOT disqualify the protein; only a TM helix in the mature-chain coordinates does).

## Stage 3 — cross-reference reconciliation

- Bd: **do reciprocal-best-hit (RBH) protein mapping as the primary path, not a fallback** — FungiDB's Bd JEL423 gene models are very likely on an older (Broad) assembly, while BFD predicted genes on `GCA_048537975.1` (V3); assume the assemblies differ until proven otherwise, and gate any coordinate-based shortcut behind an explicit assembly-version check. Flag domain calls found in both the BFD-recompute (RBH-mapped) and the curated JEL423 source as higher-confidence; calls found only in one side are reported but flagged lower-confidence.
- Bsal: same RBH mechanic but against raw AMFP13 GenBank annotation only (no curated FungiDB/InterPro source) — weaker cross-check, noted as such. Also check the GenBank annotation's own BUSCO-style **duplication** rate here, since Stage-3 is the natural place to sanity-check the Bsal protein-count anomaly from the Reference genomes section (19,449 proteins is ~2x the expected chytrid gene count) — if duplication is high, downstream Stage 4 tables must report counts per unique locus, not per transcript/isoform.

## Stage 4 — linking to metabolomics

- **Subtract media-control signal before linking.** Liquid-fraction CANOPUS classes (terpenoid/polyketide/etc.) will substantially include `*C_liq` media-blank components (Bd 1% Tryptone, Bsal 50% TGHL) — filter or background-subtract media-control features per `curated_gnps_metadata.tsv` conventions *before* any compound enters the candidate table, otherwise the linkage is measuring media chemistry, not fungal secretion.
- Source tables: `analysis/sirius_annotation/sirius_annotations.tsv` (CANOPUS/NPClassifier compound class per feature) joined to `analysis/differential_features_primary/` liquid-vs-spore contrasts, restricted to features enriched in the liquid (media) fraction, after the media-control subtraction above.
- Compound class → candidate domain family mapping: terpenoid → terpene synthase hits (PF01397/PF03936, excluding PF00494 squalene/phytoene synthase per Stage 1); polyketide → PKS hits; peptide-like/alkaloid-with-amide-bond → NRPS/RiPP-adjacent hits; unmapped classes are left out of the candidate table rather than force-fit.
- Per species/compound-class candidate table, secreted proteins (Stage 2) with matching domain evidence (Stage 1) are ranked by a **tiered/lexicographic ranking, not a weighted composite score** — a single weighted-sum score over incommensurate inputs (domain count, SignalP confidence, cluster membership, metabolite fold-change) is unfalsifiable and hides which evidence actually drove a candidate's rank. Instead: sort candidates by evidence tier (e.g. tier 1 = domain hit + cross-reference-confirmed + BGC context; tier 2 = domain hit + cross-reference-confirmed; tier 3 = domain hit only), with the metabolite's differential-abundance significance/fold-change as a tie-breaker column within a tier — and expose every evidence column in the output table so the ranking is auditable, not opaque. This is NOT real co-expression (no RNA-seq/proteomic quantification exists for Bd/Bsal life stages in this project, only a Trinity assembly used as gene-prediction evidence) — document it plainly as a heuristic ranking.

## Compute

- Stages 1–2 run inside BFD's pipeline (its own SLURM/nextflow wiring; scope with `--taxon GENUS:Batrachochytrium` so only these two genomes run) — use `-p preempt -A preempt` for any BFD job submission this work triggers, consistent with the rest of this spec.
- DeepTMHMM (standalone, this repo) requires `-p preempt_gpu -A preempt --gres=gpu:1` (GPU + singularity container, no CPU-only path).
- Stage 3 RBH mapping (BLAST/DIAMOND) and everything in Stage 4 (parsing, ranking, table-building) runs in this repo's existing pixi env via `-p preempt -A preempt` where SLURM is needed, matching the pattern already used by `analysis/*/scripts/`.

## Output

New `analysis/genome_bioactivity_linkage/` directory:
- `scripts/` — DeepTMHMM sbatch wrapper (Stage 2), RBH cross-reference script (Stage 3), and the media-subtraction + linking/ranking script (Stage 4); Stage 1–2's BFD-side tools are invoked in BFD's own tree, not duplicated here.
- Per-species domain/secretion/BGC tables pulled from BFD's outputs plus this repo's DeepTMHMM run (intermediate).
- Final linkage table(s): compound class × candidate gene × tiered evidence, one per species.
- `GENOME_BIOACTIVITY_LINKAGE.md` writeup following the `analysis/<topic>/<TOPIC>.md` convention already used by `SIRIUS_ANNOTATION.md`, including the caveats above (Bsal weaker cross-check, tiered ranking is not real co-expression, Bsal protein-count/duplication anomaly, media-control subtraction applied).

## Reproducibility / caching (SIRIUS is still incomplete)

`analysis/sirius_annotation/` is explicitly **not final** — the full native
SIRIUS run (~3,815 remaining targets) is deferred, not cancelled (see
`SIRIUS_ANNOTATION.md`), and will add annotations to `sirius_annotations.tsv`
whenever it's launched later. Stage 4 (the only stage that reads that file)
must therefore be cheap to rerun from scratch on demand, not a one-shot:

- Stages 1–3 (antiSMASH/SignalP/DeepTMHMM/PredGPI/hmmscan/cross-reference) are
  the expensive, genome-side steps and only need to run once per genome —
  each writes to a persistent, genome-keyed output directory (matching BFD's
  own persistence convention, e.g. `results/<species>/antismash_local/`,
  `results/<species>/signalp/`) and every script/sbatch step must **skip
  cleanly if its output already exists** (check-then-run, not
  overwrite-always), so re-invoking the stage-1–3 driver after an interrupted
  run or a later re-review doesn't repeat finished work.
- Stage 4 (linking) is cheap pure-Python and must be re-run in full each time
  it's invoked — no caching there, no partial-update logic. It always reads
  the current `sirius_annotations.tsv` fresh, so as soon as the full native
  SIRIUS run is merged, simply re-running the Stage 4 script picks up the new
  annotations automatically with no manual bookkeeping.
- The `GENOME_BIOACTIVITY_LINKAGE.md` writeup should note the run date/
  SIRIUS-coverage snapshot (e.g. "1,885 features annotated as of 2026-08-20
  pilot merge") each time Stage 4 is rerun, so a stale linkage table is
  never mistaken for one reflecting full SIRIUS coverage.

## Implementation approach (resolved by review)

No new Nextflow pipeline in this repo, and no hand-rolled shell reimplementation
of tools BFD already wraps. Stages 1–2 run inside BFD's existing pipeline
(`ANTISMASH_RUN`, `SIGNALP_RUN`, `RUN_PFAM`, `RUN_PREDGPI` modules, scoped
with `--taxon GENUS:Batrachochytrium`) since both genomes are already BFD
samples and the annotation belongs in BFD's shared tree. DeepTMHMM (Stage 2)
and Stages 3–4 (RBH cross-reference, linking/ranking) are plain
pixi-env Python + `sbatch` scripts in this repo's own
`analysis/genome_bioactivity_linkage/scripts/`, matching this project's
existing `analysis/<topic>/scripts/` convention — a 2-genome, ~10-task
fan-out gets no benefit from a workflow manager's parallelism/resume
features, and this repo has zero prior Nextflow usage.
