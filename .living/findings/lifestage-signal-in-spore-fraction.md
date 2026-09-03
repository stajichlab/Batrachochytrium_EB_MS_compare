## F-003: Life-stage (Zoospore-vs-Developed) signal is concentrated in the spore fraction, not the liquid
**Status:** direction supported; numbers superseded and the 2-state/3-state sub-claim RETRACTED (2026-09-02 corrected re-run)
**Claim (SUPERSEDED NUMBERS — see the 2026-09-02 amendment below for the values that stand):** On the Everything-Bagel table, collapsing Sporangium+Mature into a single Developed stage and contrasting Zoospore vs Developed within each species x matrix, the life-stage signal is concentrated in the **spore fraction**: Bd 5,638 and Bsal 7,211 FDR-significant features (n=5 vs 10) vs only 54–536 in the **liquid supernatant** (liq n=10 vs 20). The liq-vs-spore (secreted vs cellular) family is far larger still (17.6 k–27.9 k significant), dominating every life-stage effect — consistent with F-002's matrix dominance.

**Correction #1 (2026-09-02, ITSELF SUPERSEDED later the same day — its 5,507 figure comes from the unfiltered 38,547-feature table and its Bd-3-state/Bsal-2-state reading is retracted; kept for the audit trail):** the collapse was originally justified as "*every* within-matrix Sporangium-vs-Mature contrast was 0-significant". That is **false for Bd's spore fraction**: `dendrobatidis_spore_Sporangium_vs_spore_Mature` = **5,507** significant of 21,816 tested. The 0-significant statement holds only for Bsal spore (0) and for both species' liq contrasts (0, except Bd liq Zoospore-vs-Sporangium = 398). Bd's spore fraction is a genuine 3-state trajectory (Zoo-vs-Spor 3,396; Spor-vs-Mature 5,507; Zoo-vs-Mature 7,566) whereas Bsal spore really is 2-state. The concentration-in-spore-fraction claim above is unaffected, but the Bd collapse discards 5,507 significant features and should be re-derived uncollapsed or with an ordinal trend test.
**Implications:** Developmental reprogramming of the metabolome happens in/on the cells (spore fraction), while the supernatant is dominated by shared media/secreted chemistry that does not strongly re-partition across the zoospore→developed transition. For the "secreted compound / high bioactivity" goal, the informative layers are (a) liq-vs-spore (secreted_vs_cellular) contrasts and (b) life-stage hits that are also liq-enriched (`is_secreted_candidate`) — not the raw life-stage contrast in liq.
**Tags:** metabolomics, differential-abundance, life-stage, zoospore, matrix, everything-bagel, sirius, bioactivity


**AMENDED AGAIN 2026-09-02 (full corrected re-run).** Two things changed the
basis of this finding, and the headline claim survives while its numbers and
one sub-claim do not.

*What broke.* The liq groups contained 30 uninoculated `C_liq` media blanks
(50% of every liq group), so the "only 54-536 significant in liq" figure that
this finding used as evidence for a cell-associated signal was a dilution
artifact. Fungal-only, with the artifact filter and an exact permutation null,
the liq contrasts are **Bd 365** and **Bsal 2,300** significant, not 536/54 --
Bsal's went UP 43x. The asymmetry is therefore much weaker than reported.

*What holds.* The spore fraction still carries more life-stage signal
(Bd 3,219 / Bsal 3,896 vs liq 365 / 2,300), and the ordination gives
independent support: PCoA2 resolves a clean monotonic stage gradient in spore
(Zoospore +0.246 -> Sporangium -0.074 -> Mature -0.172) against a flat liq
fraction (+0.009 / -0.003 / -0.005). The direction of F-003 is intact.

*What is retracted.* The "Bsal spore is genuinely 2-state, Bd is 3-state"
reading is withdrawn entirely. At n=5v5 a BH call needs ~16% of features to
separate perfectly at once, so "0 significant" was never evidence of no
signal. Threshold-free complete-separation enrichment shows **all 12
within-matrix stage pairs carry ordered signal and not one is BH-callable**:
Bd liq 3.0x/20.3x/16.3x, Bd spore 23.8x/42.5x/24.0x, Bsal liq 7.3x/22.6x/11.8x,
Bsal spore 25.5x/35.1x/6.7x. Bd spore Sporangium-vs-Mature -- the 5,507 that
the previous amendment treated as the key datum -- reads 0 on the
artifact-filtered table and 5,507 on the unfiltered one while sitting at 24.0x
enrichment in both. Both species have a 3-state spore trajectory. Use the
ordinal trend tier, not the pairwise BH counts, for any stage claim.

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-20 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) | Bd_massspec/Batrachochytrium_MS | life-stage sig: Bd spore 5,638 / liq 536; Bsal spore 7,211 / liq 54; liq-vs-spore sig 17.6 k–27.9 k; 103,637 sig feature-rows (33,066 unique), 6,688 NPC-annotated, 3,003 secreted candidates, 1,041 bioactivity-flagged | supports |
| 2026-08-20 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) | Bd_massspec/Batrachochytrium_MS | Collapse justification: within-matrix Sporangium-vs-Mature contrasts were 0-significant (F-002 scan) → collapsing for power is empirically supported | superseded 2026-09-02 |
| 2026-09-02 | mycelium-repair-2026-09-02 | analysis/differential_features/comparison_summary.csv | Batrachochytrium_EB_MS_compare | Collapse justification is species-specific: Bd spore Sporangium-vs-Mature = 5,507 sig (not 0); Bsal spore = 0. Bd spore trajectory 3,396 / 5,507 / 7,566 across the three stage pairs | partially refutes |
| 2026-09-02 | mycelium-repair-2026-09-02 | aligned_features.csv + curated_gnps_metadata.tsv via fungal_over_blank_ratio | Batrachochytrium_EB_MS_compare | Media-blank filter never applied to the secreted/USI-curation paths: Bd 1,187/13,170 (9.0%) and Bsal 1,436/7,179 (20.0%) pass ≥2× over blank; top-100 MS² grids 2/100 and 7/100; all four named priority targets at or below blank | qualifies |

| 2026-09-02 | corrected-reanalysis-2026-09-02 | 25,157 artifact-filtered feats, 60 fungal samples, exact permutation null | Batrachochytrium_EB_MS_compare | life-stage sig: Bd spore 3,219 / liq 365; Bsal spore 3,896 / liq 2,300; all 12 stage pairs 3.0-42.5x separation-enriched, none BH-callable; trend tier Bd liq 3,114 / spore 5,027, Bsal liq 4,912 / spore 4,534 | amends + partially refutes |

### Open Questions
- Are the bioactivity-flagged hits (n=1,041, curated keyword heuristic) confirmed by direct MS/MS inspection (metabolomics-USI resolver) after controlling for PEG/adduct art facts?
- Do the spore-fraction life-stage hits (Bd 5.6 k / Bsal 7.2 k) represent conidial/spore wall remodeling vs active secretion, and what fraction are liq-enriched (secreted candidates)?
- After media-blank filtering, how much of the spore-fraction life-stage signal survives, and does the Bd 3-state / Bsal 2-state asymmetry persist?

---
