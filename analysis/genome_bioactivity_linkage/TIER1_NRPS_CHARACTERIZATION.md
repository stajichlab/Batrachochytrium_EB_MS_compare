# Tier-1 NRPS Candidate Characterization

Follow-up on `GENOME_BIOACTIVITY_LINKAGE.md`: deep characterization of the
single Tier-1 NRPS candidate gene per species (the only proteins that score
both BGC-context AND RBH cross-ref confirmed), using the antiSMASH BGC
regions and NRPSPredictor domain predictions, orthology to the NCBI
annotation, and MIBiG knownclusterblast results.

## Candidates

| Species | BFD protein | NCBI ortholog | NCBI locus | antiSMASH region | MIBiG top hit | Expressed? |
|---------|-------------|---------------|------------|------------------|---------------|------------|
| Bd JEL423 | `FCC698BD_004035-T1` | XJO74646.1 (+XJO74647.1, 2nd ORF) | BDV3_005439 | CP161927.1 rgn1 **NRPS** | hassallidin C/D (BGC0000369), anachelin (BGC0002532), gramicidin LgrC | yes (2534 rc) |
| Bsal AMFP13 | `F61BA062_001377-T1` | KAJ1345353.1 | BSLG_000866 | LYON02000001.1 rgn1 **NRPS** | nostophycin (BGC0001029), simC/cyanopeptolin (BGC0000334) | yes (739 rc) |

## Bd NRPS (BDV3_005439, CP161927.1 region 1, ~65.6 kb BGC, 29 CDS)

- antiSMASH calls region 1 as **NRPS** (protocluster NRPS, not contig-edge).
- The tier-1 protein is the cluster's large type-I NRPS: antiSMASH split the
  gene model into two ORFs (XJO74646.1 = `BDV3_005439_4100c1f6`, XJO74647.1 =
  `BDV3_005439`) that both encode the same C–A–PCP module snapshot:
  - Condensation domain (active), Adenylation domain with **NRPSPredictor
    substrate consensus = Ala** (alanine-activating), PCP/PP-binding ×2.
  - PFAM architecture (PF00550 PCP, PF00668 Condensation, PF00501
    AMP-binding, ×2 each) supports a ~2-module NRPS.
- A free-standing transaminase (`BDV3_005441`/XJO74649.1, aSDomain
  `Aminotran_1_2`; MIBiG hits to cystathionine γ-synthase/γ-lyase) sits
  immediately downstream — plausible tailoring/amino-group chemistry for the
  peptide.
- **MIBiG knownclusterblast** (all ~30% identity — distant):
  - **hassallidin C/D** BGC (BGC0000369, NRPS+saccharide antifungal
    glycolipopeptides), best hit to the cluster NRPS (E≈2e-129),
  - **anachelin** BGC (BGC0002532, NRPS+PKS siderophore) — the anachelin
    cluster also encodes a chorismate/utilization + Benzoate-CoA-ligase +
    ferrichrome-iron receptor set,
  - gramicidin S/LgrC-type synthetase (BGC0000367) and heinamide-like
    NRPS+PKS (BGC0002572).
- Interpretation: the machine is clearly a type-I NRPS (Ala-loading module),
  most similar to peptide/siderophore-type fungal/bacterial NRPS assemblies;
  nonexact homology means the exact product is **not** predictable from
  sequence alone.

## Bsal NRPS (BSLG_000866, LYON02000001.1 region 1, ~66.4 kb BGC, 15 CDS)

- antiSMASH region 1 = **NRPS**.
- Tier-1 protein = C–A–PCP module (Condensation active, Adenylation with
  substrate consensus **X** = unconfident, PCP/PP-binding); PFAM PF00501
  AMP-binding ×2 + PF00550 PCP + PF00668 Condensation.
- Co-cluster NRPS machinery is compact and split across several large ORFs:
  - `BSLG_000863` (4,383 aa; 6 PFAM domains PF17407/406/405/404/403 +
    PF03813)
  - `BSLG_000866` (target, 1,732 aa), `BSLG_000867` (3,619 aa),
    `BSLG_000868` (6,029 aa)
  - flanked by a peptidyl-prolyl cis-trans isomerase (`BSLG_000869`) and
    several hypotheticals — a small, tight NRPS island.
- **MIBiG knownclusterblast**: best hit **nostophycin** (BGC0001029,
  AEU11006.1 = NpnB, E≈7e-127, NRPS Type I) — a cyclic hexapeptide; also
  simC (BGC0000334, cyanopeptolin-like NRPS, E≈7e-59). Weaker VrtD
  (viriditoxin, iterative type-I PKS) and clavarinone cyclase hits appear
  elsewhere in the species' region-1 sets.
- Interpretation: a single-module-plus-flanking type-I NRPS most similar to
  nostophycin/cyanopeptolin cyclic-peptide assemblies; A-domain specificity
  not confidently assigned.

## Cross-species synthesis

- Both species' **only** Tier-1 biosynthetic genes are **type-I NRPS
  assemblies with siderophore/cyclopeptide-type similarity** (anachelin,
  nostophycin, cyanopeptolin are all NRP **siderophores or cyclic peptides**).
- This lines up with the metabolomics side: the dominant mapped compound
  class in both candidate tables is **"Amino acids and Peptides"** (Bd 2,574
  of 2,634 rows; Bsal 8,190 of 8,322) — i.e. the top-scoring genome hits sit
  in the same ontology family as the bulk of the liquid-enriched
  peptide-class signal. Direct candidate specifics (e.g. an Ala-signature
  hexa/heptapeptide matching an observed m/z) remain to be closed by MS².

## Outputs

- Region GBKs: `results/antismash_ncbi/{dendrobatidis/CP161927.1.region001.gbk,
  salamandrivorans/LYON02000001.1.region001.gbk}`
- MIBiG hit tables: `results/antismash_ncbi/*/knownclusterblast/region1/*_mibig_hits.html`
- Candidate tables with the new `reference_protein_id`/`ref_locus`/
  `rna_is_expressed` columns: `results/{dendrobatidis,salamandrivorans}_candidate_table.tsv`
