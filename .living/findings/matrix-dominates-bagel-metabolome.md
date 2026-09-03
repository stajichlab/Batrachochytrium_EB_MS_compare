## F-002: Matrix (liq vs spore) dominates the Everything-Bagel ordination, replicating EB
**Status:** supported (amended 2026-09-02 — blanks removed, variance shares revised down)
**Claim:** In the GNPS2 Everything-Bagel FBMN feature table (38,547 features, task e983829350de4bb39f278cbf22553247), sample separation is driven overwhelmingly by matrix (liq vs spore), not life stage or species: the all-samples PCoA axis 1 (62.9%) separates matrix, and every liq-vs-spore pairwise differential contrast within a species yields 16–24 k FDR-significant features at q<0.05 while liq-vs-liq life-stage contrasts yield ~0 (only Bd liq_Zoospore-vs-liq_Sporangium = 398). This replicates sibling EB finding F-001 (matrix dominance) on the independent Everything-Bagel feature table.
**Implications:** Any life-stage/species comparison must stratify by matrix first; the improved feature table (~9.8x more tested features than EB's 4,107-feature MZMINE3 table) preserves the same biological hierarchy — median n_tested 9.8x and median n_significant ~9.3x EB, with rank-order of comparison effect sizes essentially identical (Spearman rho 0.985 Bd / 0.996 Bsal, n=15 each).
**Tags:** metabolomics, ordination, pcoa, matrix, differential-abundance, fbmn, everything-bagel


**AMENDED 2026-09-02 (re-run on corrected data).** The original ordination
included 30 uninoculated `C_liq` media blanks among its 90 points, so its
axis 1 was partly a medium-vs-pellet axis rather than a biological one. Re-run
on the 60 fungal samples only, with the artifact filter applied (25,157
features), **the conclusion holds and is if anything cleaner**: matrix
*completely* separates on PCoA1 with no overlap (liq range [0.207, 0.320],
spore [-0.370, -0.126]), and both species behave identically (Bd liq +0.303 /
spore -0.278; Bsal +0.271 / -0.297). Variance explained is lower than
originally reported because the blanks were inflating it: axis1 **57.5%**
(was 62.9%) all-samples, **68.6%** Bd (was 75.9%), **64.8%** Bsal (was 70.6%).

PCoA2 additionally resolves a clean, monotonic life-stage gradient **in the
spore fraction only** — Zoospore +0.246, Sporangium -0.074, Mature -0.172 —
against a flat liq fraction (+0.009 / -0.003 / -0.005). That is independent
ordination-level support for F-003, and it is visible only after the blanks
are removed.

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-09-02 | corrected-reanalysis-2026-09-02 | gnps2-everything-bagel (25,157 artifact-filtered feats, 60 fungal samples) | Batrachochytrium_EB_MS_compare | PCoA1 57.5% all / 68.6% Bd / 64.8% Bsal; matrix COMPLETELY separates on axis1 (no overlap); PCoA2 monotonic spore stage gradient, flat in liq | supports (amended) |
| 2026-08-19 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) + EB (4,107 feats) | Bd_massspec/Batrachochytrium_MS | PCoA axis1 62.9% (all), 75.9%/70.6% (per-species); liq-vs-spore sig 16–24k, liq-vs-liq ~0; Spearman rho vs EB 0.985/0.996 | supports (replicates EB) |

### Open Questions
- What fraction of the 9.8x feature increase is real biology vs Everything-Bagel alignment artifacts (isobaric collapse, isotopes, adducts)? Requires crossing to the network/NPLC edges + SIRIUS annotations.
- Do the top-ranking differential features (by q-value) overlap in identity (m/z+RT+MS2-cosine) with EB's top features, or are they a new set entirely?

---
