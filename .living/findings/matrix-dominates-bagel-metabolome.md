## F-002: Matrix (liq vs spore) dominates the Everything-Bagel ordination, replicating EB
**Status:** supported
**Claim:** In the GNPS2 Everything-Bagel FBMN feature table (38,547 features, task e983829350de4bb39f278cbf22553247), sample separation is driven overwhelmingly by matrix (liq vs spore), not life stage or species: the all-samples PCoA axis 1 (62.9%) separates matrix, and every liq-vs-spore pairwise differential contrast within a species yields 16–24 k FDR-significant features at q<0.05 while liq-vs-liq life-stage contrasts yield ~0 (only Bd liq_Zoospore-vs-liq_Sporangium = 398). This replicates sibling EB finding F-001 (matrix dominance) on the independent Everything-Bagel feature table.
**Implications:** Any life-stage/species comparison must stratify by matrix first; the improved feature table (~9.8x more tested features than EB's 4,107-feature MZMINE3 table) preserves the same biological hierarchy — median n_tested 9.8x and median n_significant ~9.3x EB, with rank-order of comparison effect sizes essentially identical (Spearman rho 0.985 Bd / 0.996 Bsal, n=15 each).
**Tags:** metabolomics, ordination, pcoa, matrix, differential-abundance, fbmn, everything-bagel

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-19 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) + EB (4,107 feats) | Bd_massspec/Batrachochytrium_MS | PCoA axis1 62.9% (all), 75.9%/70.6% (per-species); liq-vs-spore sig 16–24k, liq-vs-liq ~0; Spearman rho vs EB 0.985/0.996 | supports (replicates EB) |

### Open Questions
- What fraction of the 9.8x feature increase is real biology vs Everything-Bagel alignment artifacts (isobaric collapse, isotopes, adducts)? Requires crossing to the network/NPLC edges + SIRIUS annotations.
- Do the top-ranking differential features (by q-value) overlap in identity (m/z+RT+MS2-cosine) with EB's top features, or are they a new set entirely?

---
