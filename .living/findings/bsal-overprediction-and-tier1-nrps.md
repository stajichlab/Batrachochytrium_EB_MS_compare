## F-004: Both species' top biosynthetic gene set is NRPS — siderophore/cyclopeptide-type — matching the dominant peptide mass-spec class, and Bsal's BFD proteome is substantially over-predicted
**Status:** supported
**Claim:** (a) The single Tier-1 NRPS candidate gene in each species (Bd BDV3_005439 / Bsal BSLG_000866) is a type-I NRPS whose module snapshot (C–A–PCP, Ala vs unknown substrate) and antiSMASH-knownclusterblast hits (Bd: hassallidin C/D BGC0000369, anachelin BGC0002532; Bsal: nostophycin BGC0001029, cyanopeptolin-like BGC0000334) place it among peptide/siderophore-type NRPS assemblies — at ~30% identity (distant), so product identity is not determined. This is ontologically concordant with the dominant mapped compound class in both candidate tables ("Amino acids and Peptides"). (b) Bsal's BFD proteome is over-predicted: 19,449 BFD loci vs 10,867 NCBI genes (≥2.3× over Bd's ratio), one transcript per locus (no isoform inflation), all-vs-all self-BLAST shows 30.7% near-duplicate proteins (Bd 14.7%) plus a heavier short-protein tail (Bsal 29.5% <200 aa vs Bd 18.4%); the DeepTMHMM-crashing 4,777-aa outlier is a real single-exon hypothetical ORF in BOTH BFD and NCBI annotations (BSLG_008696/KAJ1332392.1) with no domains/repeats and is transcribed.
**Implications:** Follow-up effort should concentrate on the two NRPS BGCs (both transcribed per public RNA-seq) and their co-cluster tailoring genes, validating candidate peptides against the top liq-enriched SIRIUS-annotated peptide features by MS² (USI scaffold built). Bsal candidate counts should be read knowing ~31% of its 19,449 predicted proteins may be over-predicted duplicates/short ORFs.
**Tags:** genome-bioactivity-linkage, nrps, bgc, duplication, over-prediction, sirius, peptide, metabolomics, everything-bagel

### Evidence Ledger
| Date | Run/Session | Dataset | Project | Result | Direction |
|------|-------------|---------|---------|--------|-----------|
| 2026-08-26 | gbl-followup-2026-08-26 | antiSMASH rgn1 GBK + knownclusterblast MIBiG + RBH | Bd_massspec/Batrachochytrium_MS | Bd tier-1 NRPS = C-A(Ala)-PCP + transaminase, MIBiG hassallidin/anachelin-type; Bsal = C-A(X)-PCP + flanking NRPS ORFs, MIBiG nostophycin/cyanopeptolin-type | supports |
| 2026-08-26 | gbl-followup-2026-08-26 | DIAMOND self-blastp (results/duplication_check) | Bd_massspec/Batrachochytrium_MS | Bd 8,396 proteins → 11.5% strict / 14.7% near-dups; Bsal 19,449 → 17.8% strict / 30.7% near-dups; BFD 19,449 loci vs NCBI 10,867 genes, single isoform/locus | supports |
| 2026-08-26 | gbl-followup-2026-08-26 | featureCounts + gffread GTF + run_duplication_check.sh | Bd_massspec/Batrachochytrium_MS | Outlier F61BA062_016014-T1 present in BFD and NCBI (BSLG_008696, 4,681 aa NCBI), no PFAM, no repeats, transcribed (274/179 counts) | supports |

### Open Questions
- Does either chytrid's peptide profile (top liq-enriched "Amino acids and Peptides" features) match a predicted Ala-rich / cyclo-peptide product from these NRPS clusters (SIRIUS + USI-level verification)?
- Is Bsal's over-prediction concentrated in specific gene families (segmental duplicates?) that would change the biological interpretation of its larger candidate table?

---

