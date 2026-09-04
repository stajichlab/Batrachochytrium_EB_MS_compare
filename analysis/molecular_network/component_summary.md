# Molecular-network components

- network: **9600 edges, 4423 nodes, 659 components** (largest 97)
- after the artifact filter and analysis-matrix restriction: **599 components**, largest **61**, 275 with >= 3 members
- background blank-clearing rate: dendrobatidis 4.39%, salamandrivorans 7.40%

## dendrobatidis

- components significantly blank-enriched (q<0.05, >=3 members): **11**
- components ENTIRELY blank-clearing (>=3 members): **1**

| component | size | blank-clearing | q | dominant class | examples |
|---|---|---|---|---|---|
| 0 | 59 | 25 | 3.92e-16 | Fatty acids | [2-Acetyloxy-4-(3-hexadecoxy-2-hydroxypropoxy)-4-oxobutyl]-trimethylaz |
| 209 | 57 | 9 | 4.31e-02 | Fatty acids | (2-aminoethoxy)[2-hydroxy-3-(octadeca-9,12,15-trienoyloxy)propoxy]phos |
| 369 | 29 | 11 | 3.72e-06 | Fatty acids | Chrysogeside A; (Z)-N-[(E)-3-hydroxy-1-[3,4,5-trihydroxy-6-(hydroxymet |
| 95 | 17 | 7 | 4.89e-04 | Amino acids and Peptides | 6,9-dibenzyl-12-(butan-2-yl)-3,15-bis(hydroxymethyl)-1,4,7,10,13,16-he |
| 408 | 14 | 10 | 6.56e-09 | Fatty acids | 2-{[2-({6-[(2-dodecanamido-3-hydroxyoctadec-4-en-1-yl)oxy]-4,5-dihydro |
| 52 | 9 | 4 | 2.92e-02 | Amino acids and Peptides | 2-[[1-[3-Phenyl-2-(pyrrolidine-2-carbonylamino)propanoyl]pyrrolidine-2 |
| 381 | 6 | 4 | 5.16e-03 | Amino acids and Peptides | 6-amino-2-{2-[2-(2-amino-3-carbamoylpropanamido)-4-carboxybutanamido]p |
| 118 | 5 | 3 | 4.31e-02 | Amino acids and Peptides | tert-butyl N-[(2S)-1-[(2S)-2-[[4-[[2-(benzylamino)-2-oxoethyl]amino]-1 |
| 361 | 5 | 3 | 4.31e-02 | Fatty acids | [(2R)-3-[2-aminoethoxy(hydroxy)phosphoryl]oxy-2-nonanoyloxypropyl] (E) |
| 339 | 4 | 4 | 4.89e-04 | Fatty acids | [5-[2-[2-[2-(3-Methoxyphenyl)ethyl]phenoxy]ethyl]-1-methylpyrrolidin-3 |

## salamandrivorans

- components significantly blank-enriched (q<0.05, >=3 members): **12**
- components ENTIRELY blank-clearing (>=3 members): **2**

| component | size | blank-clearing | q | dominant class | examples |
|---|---|---|---|---|---|
| 0 | 59 | 15 | 3.67e-03 | Fatty acids | [2-Acetyloxy-4-(3-hexadecoxy-2-hydroxypropoxy)-4-oxobutyl]-trimethylaz |
| 154 | 35 | 19 | 2.33e-10 | Amino acids and Peptides | H-DL-Leu-Gly-DL-xiThr-DL-xiIle-DL-Pro-Gly-OH; H-Val-Val-Pro-Pro-Phe-OH |
| 369 | 29 | 10 | 3.83e-03 | Fatty acids | Chrysogeside A; (Z)-N-[(E)-3-hydroxy-1-[3,4,5-trihydroxy-6-(hydroxymet |
| 338 | 19 | 8 | 3.83e-03 | Fatty acids | 9-[1-(2-Aminoethoxy)-3-heptadecanoyloxy-1,2-dihydroxypropan-2-yl]oxy-9 |
| 95 | 17 | 6 | 4.94e-02 | Amino acids and Peptides | 6,9-dibenzyl-12-(butan-2-yl)-3,15-bis(hydroxymethyl)-1,4,7,10,13,16-he |
| 346 | 12 | 6 | 8.76e-03 | Amino acids and Peptides | H-Pro-Ser-Pro-Ser-Pro-Ser-al; 2-[2-(2-amino-3-phenylpropanamido)-3-hyd |
| 52 | 9 | 5 | 1.62e-02 | Amino acids and Peptides | 2-[[1-[3-Phenyl-2-(pyrrolidine-2-carbonylamino)propanoyl]pyrrolidine-2 |
| 417 | 8 | 5 | 8.76e-03 | Amino acids and Peptides | Citrusin Iii; Citrusin Iii; H-DL-xiIle-DL-Leu-DL-Val-DL-Tyr-OH |
| 42 | 7 | 6 | 3.21e-04 | Amino acids and Peptides | H-Ala-Pro-Glu-Ala-Val-OH; Gln-Glu-Pro-Val-Leu; 15-benzyl-21-(butan-2-y |
| 197 | 6 | 4 | 2.21e-02 | Amino acids and Peptides | H-Pro-Ser-Pro-Ser-Pro-Ser-al; H-Ala-Pro-Asp-Cys-Arg-Pro-al; L-Iditol,  |

## DeltaMZ homologous steps (both ends surviving the artifact filter)

| step | edges | median cosine |
|---|---|---|
| 2xCH2 | 156 | 0.914 |
| CH2 (homolog) | 152 | 0.905 |
| H2 (saturation) | 126 | 0.876 |
| O (oxidation) | 97 | 0.886 |
| H2O | 44 | 0.863 |
| C2H2 | 40 | 0.862 |
| NH | 4 | 0.806 |

Total classified edges: 619 of 4703 intra-matrix edges (13.2%).
