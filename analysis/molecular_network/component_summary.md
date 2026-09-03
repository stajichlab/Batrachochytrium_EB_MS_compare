# Molecular-network components

- network: **9600 edges, 4423 nodes, 659 components** (largest 97)
- after the artifact filter and analysis-matrix restriction: **551 components**, largest **47**, 226 with >= 3 members
- background blank-clearing rate: dendrobatidis 4.16%, salamandrivorans 7.28%

## dendrobatidis

- components significantly blank-enriched (q<0.05, >=3 members): **5**
- components ENTIRELY blank-clearing (>=3 members): **0**

| component | size | blank-clearing | q | dominant class | examples |
|---|---|---|---|---|---|
| 0 | 31 | 8 | 8.10e-03 | Fatty acids | 4-[3-hexadecanoyloxy-2-[(Z)-octadec-9-enoyl]oxypropoxy]-2-(trimethylaz |
| 369 | 19 | 6 | 1.60e-02 | Fatty acids | (Z)-N-[(E)-3-hydroxy-1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]ox |
| 95 | 7 | 6 | 1.91e-05 | Amino acids and Peptides | 2-amino-9-[(3R,4S,5R)-3,4-dihydroxy-5-(hydroxymethyl)oxolan-2-yl]-2,3- |
| 471 | 4 | 3 | 3.07e-02 | — | — |
| 361 | 4 | 3 | 3.07e-02 | Fatty acids | [(2R)-3-[2-aminoethoxy(hydroxy)phosphoryl]oxy-2-nonanoyloxypropyl] (E) |

## salamandrivorans

- components significantly blank-enriched (q<0.05, >=3 members): **8**
- components ENTIRELY blank-clearing (>=3 members): **2**

| component | size | blank-clearing | q | dominant class | examples |
|---|---|---|---|---|---|
| 154 | 26 | 11 | 2.68e-04 | Amino acids and Peptides | H-DL-Leu-Gly-DL-xiThr-DL-xiIle-DL-Pro-Gly-OH; H-Val-Val-Pro-Pro-Phe-OH |
| 369 | 19 | 7 | 2.66e-02 | Fatty acids | (Z)-N-[(E)-3-hydroxy-1-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]ox |
| 42 | 7 | 6 | 2.68e-04 | Amino acids and Peptides | H-Ala-Pro-Glu-Ala-Val-OH; Gln-Glu-Pro-Val-Leu; 15-benzyl-21-(butan-2-y |
| 417 | 6 | 4 | 2.66e-02 | Amino acids and Peptides | Citrusin Iii; Citrusin Iii; H-DL-xiIle-DL-Leu-DL-Val-DL-Tyr-OH |
| 52 | 6 | 4 | 2.66e-02 | Amino acids and Peptides | 2-[[1-[3-Phenyl-2-(pyrrolidine-2-carbonylamino)propanoyl]pyrrolidine-2 |
| 197 | 6 | 4 | 2.66e-02 | Amino acids and Peptides | H-Pro-Ser-Pro-Ser-Pro-Ser-al; H-Ala-Pro-Asp-Cys-Arg-Pro-al; L-Iditol,  |
| 468 | 3 | 3 | 2.66e-02 | Amino acids and Peptides | 12-benzyl-3,15-bis(1-hydroxyethyl)-6-[(4-hydroxyphenyl)methyl]-9-(2-me |
| 646 | 3 | 3 | 2.66e-02 | Terpenoids | N''-(2-{28'-amino-4-[8-(ethylamino)-8-hydroxy-11-methyl-7-methylidene- |

## DeltaMZ homologous steps (both ends surviving the artifact filter)

| step | edges | median cosine |
|---|---|---|
| CH2 (homolog) | 68 | 0.884 |
| 2xCH2 | 49 | 0.885 |
| H2 (saturation) | 49 | 0.877 |
| O (oxidation) | 44 | 0.887 |
| C2H2 | 13 | 0.842 |
| H2O | 11 | 0.863 |
| NH | 3 | 0.866 |

Total classified edges: 237 of 2781 intra-matrix edges (8.5%).
