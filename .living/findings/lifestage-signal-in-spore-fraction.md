## F-003: Life-stage (Zoospore-vs-Developed) signal is concentrated in the spore fraction, not the liquid
**Status:** supported
**Claim:** On the Everything-Bagel table, collapsing Sporangium+Mature into a single Developed stage (justified by F-002's ~0 within-matrix stage contrasts) and contrasting Zoospore vs Developed within each species x matrix, the life-stage signal is concentrated in the **spore fraction**: Bd 5,638 and Bsal 7,211 FDR-significant features (n=5 vs 10) vs only 54–536 in the **liquid supernatant** (liq n=10 vs 20). The liq-vs-spore (secreted vs cellular) family is far larger still (17.6 k–27.9 k significant), dominating every life-stage effect — consistent with F-002's matrix dominance.
**Implications:** Developmental reprogramming of the metabolome happens in/on the cells (spore fraction), while the supernatant is dominated by shared media/secreted chemistry that does not strongly re-partition across the zoospore→developed transition. For the "secreted compound / high bioactivity" goal, the informative layers are (a) liq-vs-spore (secreted_vs_cellular) contrasts and (b) life-stage hits that are also liq-enriched (`is_secreted_candidate`) — not the raw life-stage contrast in liq.
**Tags:** metabolomics, differential-abundance, life-stage, zoospore, matrix, everything-bagel, sirius, bioactivity

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-20 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) | Bd_massspec/Batrachochytrium_MS | life-stage sig: Bd spore 5,638 / liq 536; Bsal spore 7,211 / liq 54; liq-vs-spore sig 17.6 k–27.9 k; 103,637 sig feature-rows (33,066 unique), 6,688 NPC-annotated, 3,003 secreted candidates, 1,041 bioactivity-flagged | supports |
| 2026-08-20 | bagel-pipeline-2026-08-19 | gnps2-everything-bagel (38,547 feats) | Bd_massspec/Batrachochytrium_MS | Collapse justification: within-matrix Sporangium-vs-Mature contrasts were 0-significant (F-002 scan) → collapsing for power is empirically supported | supports |

### Open Questions
- Are the bioactivity-flagged hits (n=1,041, curated keyword heuristic) confirmed by direct MS/MS inspection (metabolomics-USI resolver) after controlling for PEG/adduct art facts?
- Do the spore-fraction life-stage hits (Bd 5.6 k / Bsal 7.2 k) represent conidial/spore wall remodeling vs active secretion, and what fraction are liq-enriched (secreted candidates)?

---
