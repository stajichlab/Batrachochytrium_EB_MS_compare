# Session Log Registry

| Date | Session ID | Project | Branch | Duration | Files Changed | Summary | Key Outputs | Status | Tags | Log |
|------|-----------|---------|--------|----------|---------------|---------|-------------|--------|------|-----|
| 2026-08-19 | bagel-pipeline-2026-08-19 | Bd_massspec/Batrachochytrium_MS | main | ~40 min | pixi.toml, analysis/ordination/scripts/build_ordination_table.py, figures, differential_features/, DATA_MANIFEST.md, .living/{findings,learnings,decisions} | Set up pixi env; ported EB ordination+differential pipeline to the Everything-Bagel feature table (38,547 feats); ran all 3 scripts; compared vs EB (Spearman rho 0.985/0.996, ~9.8x feature scale-up) | figures/pcoa_all* + by_species/; 30 differential dirs w/ volcano+top_features; comparison_summary.csv; F-002 finding | complete | porting, ordination, differential | |
