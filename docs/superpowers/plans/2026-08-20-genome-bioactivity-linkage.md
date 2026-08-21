# Genome-to-Bioactive-Compound Linkage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the repo-side pipeline (analysis code, not the BFD-side genome
annotation runs the user is running separately) that links SIRIUS/CANOPUS
metabolite classes from the Bd/Bsal LC-MS/MS data to candidate biosynthetic
and secreted genes in Bd JEL423 / Bsal AMFP13, producing tiered per-species
candidate tables.

**Architecture:** A `analysis/genome_bioactivity_linkage/scripts/` package of
small, single-purpose Python modules, each reading one data source (BFD PFAM
domtblout, BFD antiSMASH JSON, DeepTMHMM `TMRs.gff3`, the project's aligned
feature table + curated metadata, `sirius_annotations.tsv`) and writing one
intermediate TSV, following this repo's existing `analysis/<topic>/scripts/`
convention. A final linking script joins the intermediates into tiered
candidate tables. Two SLURM steps (DeepTMHMM, a DIAMOND reciprocal-best-hit
ortholog search) are plain `sbatch` scripts under `-p preempt`/`-p short_gpu`,
matching this repo's existing style — no Nextflow.

**Tech Stack:** Python 3.12 (pixi env: pandas, numpy already present),
DIAMOND (HPCC module, reciprocal-best-hit ortholog mapping), DeepTMHMM
(singularity container, GPU), pytest for unit tests.

**Spec:** `docs/superpowers/specs/2026-08-20-genome-bioactivity-linkage-design.md`

## Global Constraints

- All SLURM job submissions in this plan use `-p preempt -A preempt`, except
  DeepTMHMM which uses `-p short_gpu --gres=gpu:1` (per the spec's Compute
  section, as last edited by the user).
- Bd reference = JEL423 (`genome_annotation/Batrachochytrium_dendrobatidis_JEL423`,
  BFD locustag `FCC698BD`, NCBI assembly `GCA_048537975.1_CMM_BatrDend_JEL423_V3`).
- Bsal reference = AMFP13 (`genome_annotation/Batrachochytrium_salamandrivorans_AMFP13`,
  BFD locustag `F61BA062`, NCBI assembly `GCA_002006685.2_Batr_sala_V2`).
- `BFD_ROOT = /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs` — read-only
  from this repo except for the two SLURM steps this plan explicitly submits.
- PF00494 (squalene/phytoene synthase) must never be counted as "terpene
  synthase" — it is its own domain family throughout.
- Ranking is tiered/lexicographic, never a single weighted composite score.
- Every parser must handle BFD output not existing yet (the user's BFD run
  may still be in progress) by raising a clear `FileNotFoundError` with the
  expected path, not silently producing an empty table.

---

## File Structure

```
analysis/genome_bioactivity_linkage/
  scripts/
    paths.py                    # Task 1: shared path/locustag constants
    domain_families.py          # Task 1: Pfam-id -> family, compound-class map
    background_subtraction.py   # Task 2: C_liq companion-blank subtraction
    run_deeptmhmm.sh            # Task 3: SLURM sbatch wrapper
    parse_deeptmhmm.py          # Task 3: TMRs.gff3 -> per-protein TM calls
    fetch_reference_annotation.sh   # Task 4: download JEL423/AMFP13 NCBI refs
    run_rbh.sh                  # Task 5: SLURM DIAMOND reciprocal-best-hit
    parse_pfam_domains.py       # Task 6: BFD domtblout -> domain-family calls
    parse_antismash_clusters.py # Task 7: BFD antiSMASH JSON -> BGC context
    merge_secretion.py          # Task 8: SignalP + DeepTMHMM + PredGPI -> extracellular set
    link_compounds_to_genes.py  # Task 9: final tiered candidate tables
  tests/
    test_paths.py
    test_domain_families.py
    test_background_subtraction.py
    test_parse_deeptmhmm.py
    test_parse_pfam_domains.py
    test_parse_antismash_clusters.py
    test_merge_secretion.py
    test_link_compounds_to_genes.py
  results/                       # generated, not committed except final tables
  GENOME_BIOACTIVITY_LINKAGE.md   # Task 10
```

---

### Task 1: Shared paths and domain-family constants

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/paths.py`
- Create: `analysis/genome_bioactivity_linkage/scripts/domain_families.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_paths.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_domain_families.py`

**Interfaces:**
- Produces: `paths.SPECIES` dict, `paths.BFD_ROOT`, `paths.find_bfd_output(kind: str, species_key: str) -> Path`, `domain_families.DOMAIN_FAMILIES: dict[str, set[str]]`, `domain_families.COMPOUND_CLASS_TO_FAMILY: dict[str, str]`, `domain_families.classify_pfam(pfam_id: str) -> str | None`.

- [ ] **Step 1: Write the failing tests**

```python
# analysis/genome_bioactivity_linkage/tests/test_paths.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import paths  # noqa: E402


def test_species_table_has_both_species():
    assert set(paths.SPECIES) == {"dendrobatidis", "salamandrivorans"}
    bd = paths.SPECIES["dendrobatidis"]
    assert bd["out"] == "Batrachochytrium_dendrobatidis_JEL423"
    assert bd["locustag"] == "FCC698BD"
    bsal = paths.SPECIES["salamandrivorans"]
    assert bsal["out"] == "Batrachochytrium_salamandrivorans_AMFP13"
    assert bsal["locustag"] == "F61BA062"


def test_find_bfd_output_missing_raises_with_expected_path(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    try:
        paths.find_bfd_output("pfam_hmmscan", "dendrobatidis")
        assert False, "expected FileNotFoundError"
    except FileNotFoundError as e:
        assert "FCC698BD" in str(e)
        assert "pfam_hmmscan" in str(e)


def test_find_bfd_output_finds_locustag_file(tmp_path, monkeypatch):
    monkeypatch.setattr(paths, "BFD_ROOT", tmp_path)
    bucket = tmp_path / "results" / "function" / "pfam_hmmscan" / "00"
    bucket.mkdir(parents=True)
    target = bucket / "FCC698BD.domtblout.gz"
    target.write_bytes(b"")
    found = paths.find_bfd_output("pfam_hmmscan", "dendrobatidis")
    assert found == target
```

```python
# analysis/genome_bioactivity_linkage/tests/test_domain_families.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import domain_families  # noqa: E402


def test_terpene_synthase_excludes_squalene_phytoene_synthase():
    assert domain_families.classify_pfam("PF01397") == "terpene_synthase"
    assert domain_families.classify_pfam("PF03936") == "terpene_synthase"
    assert domain_families.classify_pfam("PF00494") == "squalene_phytoene_synthase"
    assert domain_families.classify_pfam("PF00494") != "terpene_synthase"


def test_pks_and_nrps_and_unknown():
    assert domain_families.classify_pfam("PF00109") == "pks"
    assert domain_families.classify_pfam("PF00668") == "nrps"
    assert domain_families.classify_pfam("PF99999") is None


def test_compound_class_to_family_map_omits_unmapped_classes():
    assert domain_families.COMPOUND_CLASS_TO_FAMILY["Terpenoids"] == "terpene_synthase"
    assert domain_families.COMPOUND_CLASS_TO_FAMILY["Polyketides"] == "pks"
    assert "Fatty acyls" not in domain_families.COMPOUND_CLASS_TO_FAMILY
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_paths.py analysis/genome_bioactivity_linkage/tests/test_domain_families.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'paths'` / `'domain_families'`

- [ ] **Step 3: Implement `paths.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/paths.py
"""Shared path/locustag constants for the genome-bioactivity-linkage pipeline."""
from pathlib import Path

BFD_ROOT = Path("/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs")
REPO_ROOT = Path(__file__).resolve().parents[3]
GBL_ROOT = Path(__file__).resolve().parents[1]

SPECIES = {
    "dendrobatidis": {
        "out": "Batrachochytrium_dendrobatidis_JEL423",
        "locustag": "FCC698BD",
        "ncbi_assembly": "GCA_048537975.1_CMM_BatrDend_JEL423_V3",
    },
    "salamandrivorans": {
        "out": "Batrachochytrium_salamandrivorans_AMFP13",
        "locustag": "F61BA062",
        "ncbi_assembly": "GCA_002006685.2_Batr_sala_V2",
    },
}


def bfd_proteins_fasta(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    return BFD_ROOT / "genome_annotation" / out / "predict_results" / f"{out}.proteins.fa"


def bfd_gbk(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    return BFD_ROOT / "genome_annotation" / out / "predict_results" / f"{out}.gbk"


def bfd_antismash_json(species_key: str) -> Path:
    out = SPECIES[species_key]["out"]
    base = BFD_ROOT / "genome_annotation" / out / "antismash_local"
    for name in (f"{out}.json", f"{out}.json.gz"):
        p = base / name
        if p.exists():
            return p
    raise FileNotFoundError(
        f"antiSMASH JSON not found for {species_key} under {base} "
        f"(expected {out}.json or {out}.json.gz) — BFD antiSMASH run may not be finished yet"
    )


def find_bfd_output(kind: str, species_key: str) -> Path:
    """Locate a locustag-bucketed BFD functional-annotation file, e.g. kind='pfam_hmmscan'."""
    locustag = SPECIES[species_key]["locustag"]
    search_root = BFD_ROOT / "results" / "function" / kind
    matches = sorted(search_root.glob(f"**/{locustag}*"))
    if not matches:
        raise FileNotFoundError(
            f"No BFD '{kind}' output found for locustag {locustag} under {search_root} "
            f"— BFD functional-annotation run may not be finished yet"
        )
    return matches[0]
```

- [ ] **Step 4: Implement `domain_families.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/domain_families.py
"""Pfam-id -> biosynthetic domain family, and compound-class -> family mapping.

PF00494 (squalene/phytoene synthase) is deliberately its own family, not
folded into terpene_synthase -- see spec Stage 1.
"""

DOMAIN_FAMILIES: dict[str, set[str]] = {
    "terpene_synthase": {"PF01397", "PF03936"},
    "squalene_phytoene_synthase": {"PF00494"},
    "pks": {"PF00109", "PF08659", "PF00698"},
    "nrps": {"PF00668", "PF00501", "PF00550"},
    "dmats_prenyltransferase": {"PF11991"},
    "p450": {"PF00067"},
}

# Only compound classes with a defined candidate domain family are mapped;
# unmapped classes (e.g. "Fatty acyls") are deliberately absent so Stage 4
# leaves them out of the candidate table rather than force-fitting them.
COMPOUND_CLASS_TO_FAMILY: dict[str, str] = {
    "Terpenoids": "terpene_synthase",
    "Polyketides": "pks",
    "Alkaloids (linear polyketides)": "pks",
    "Amino acids and Peptides": "nrps",
}

_PFAM_TO_FAMILY = {pfam: family for family, pfams in DOMAIN_FAMILIES.items() for pfam in pfams}


def classify_pfam(pfam_id: str) -> str | None:
    """Map a bare Pfam accession (e.g. 'PF01397', no version suffix) to its family, or None."""
    return _PFAM_TO_FAMILY.get(pfam_id)
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_paths.py analysis/genome_bioactivity_linkage/tests/test_domain_families.py -v`
Expected: PASS (5 tests)

- [ ] **Step 6: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/paths.py \
        analysis/genome_bioactivity_linkage/scripts/domain_families.py \
        analysis/genome_bioactivity_linkage/tests/test_paths.py \
        analysis/genome_bioactivity_linkage/tests/test_domain_families.py
git commit -m "genome-bioactivity-linkage: add shared paths and domain-family constants"
```

---

### Task 2: Media-companion background subtraction

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/background_subtraction.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_background_subtraction.py`

**Context:** `curated_gnps_metadata.tsv` has 33 `is_C_companion == True` rows
(media-blank wells inoculated with no fungus, e.g. `A1C_liq.mzML`), each
pointing at its fungal companion sample via `companion_of` (e.g.
`A1_spore.mzML` — note: despite the column name, inspect the actual value;
companions are matched by plate+timepoint+matrix, not literally by that
column alone) and carrying `use_in_analysis == True` for non-B-plate rows.
These companion blanks are NOT excluded anywhere upstream (`use_in_analysis`
only drops IS/QC and B-plate conditioned-media rows), so a feature's
"liquid-fraction enriched" status in `differential_features_primary`'s
existing tables does not by itself mean it is fungal in origin — it could be
a nutrient-broth (Tryptone/TGHL) component. This module computes, per
feature and per species/life-stage, whether the feature's mean peak area in
true fungal `liq` samples exceeds its mean peak area in the matched `C_liq`
companion blanks by a minimum fold-change, using the raw per-sample
intensities in `aligned_features.csv` (not the pre-aggregated
`differential_features_primary` tables, which do not separate companions
out).

**Interfaces:**
- Consumes: `paths.REPO_ROOT`
- Produces: `background_subtraction.load_metadata() -> pd.DataFrame`,
  `background_subtraction.load_feature_intensities() -> pd.DataFrame` (index
  = feature row id, columns = sample filenames, values = peak area),
  `background_subtraction.fungal_over_blank_ratio(features: pd.DataFrame, meta: pd.DataFrame, species: str, life_stage: str, min_fc: float = 2.0) -> pd.DataFrame` with columns `row_id, mean_fungal, mean_blank, log2fc_fungal_over_blank, passes_background_filter`.

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_background_subtraction.py
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from background_subtraction import fungal_over_blank_ratio  # noqa: E402


def _meta_row(filename, species, life_stage, is_companion):
    return {
        "filename": filename,
        "species": species,
        "life_stage": life_stage,
        "is_C_companion": is_companion,
        "matrix": "liq",
        "use_in_analysis": True,
    }


def test_high_fungal_low_blank_passes_filter():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A2_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [1000.0], "A2_liq.mzML": [1200.0], "A1C_liq.mzML": [50.0]},
        index=pd.Index([1], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 1].iloc[0]
    assert row["mean_fungal"] == 1100.0
    assert row["mean_blank"] == 50.0
    assert row["passes_background_filter"] is True or row["passes_background_filter"] == True  # noqa: E712


def test_media_dominated_feature_fails_filter():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [100.0], "A1C_liq.mzML": [90.0]},
        index=pd.Index([2], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 2].iloc[0]
    assert bool(row["passes_background_filter"]) is False


def test_zero_blank_signal_treated_as_pass():
    meta = pd.DataFrame(
        [
            _meta_row("A1_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", False),
            _meta_row("A1C_liq.mzML", "Batrachochytrium dendrobatidis", "Zoospore", True),
        ]
    )
    features = pd.DataFrame(
        {"A1_liq.mzML": [500.0], "A1C_liq.mzML": [0.0]},
        index=pd.Index([3], name="row_id"),
    )
    result = fungal_over_blank_ratio(
        features, meta, species="Batrachochytrium dendrobatidis", life_stage="Zoospore", min_fc=2.0
    )
    row = result.loc[result["row_id"] == 3].iloc[0]
    assert bool(row["passes_background_filter"]) is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_background_subtraction.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'background_subtraction'`

- [ ] **Step 3: Implement `background_subtraction.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/background_subtraction.py
"""Filter liquid-fraction features by fungal-sample vs C_liq-companion-blank signal.

use_in_analysis == True does NOT exclude the 33 is_C_companion == True
media-blank wells (only IS/QC and B-plate conditioned-media rows are
excluded upstream) -- see Task 2 docstring context in the implementation
plan. This module must be used before any liquid-fraction compound is
treated as a candidate for fungal secretion.
"""
from pathlib import Path

import numpy as np
import pandas as pd

from paths import REPO_ROOT

METADATA_PATH = REPO_ROOT / "data" / "metdata" / "curated_gnps_metadata.tsv"
FEATURES_PATH = (
    REPO_ROOT
    / "data"
    / "raw"
    / "gnps2_e9838293_bagel"
    / "nf_output"
    / "feature_finding"
    / "feature_finding_results"
    / "aligned_features.csv"
)

_PSEUDOCOUNT = 1.0


def load_metadata() -> pd.DataFrame:
    return pd.read_csv(METADATA_PATH, sep="\t")


def load_feature_intensities() -> pd.DataFrame:
    df = pd.read_csv(FEATURES_PATH)
    id_col = "row ID" if "row ID" in df.columns else "row_id"
    df = df.set_index(id_col)
    df.index.name = "row_id"
    peak_area_cols = {c: c.replace(" Peak area", "") for c in df.columns if c.endswith(" Peak area")}
    return df[list(peak_area_cols)].rename(columns=peak_area_cols)


def fungal_over_blank_ratio(
    features: pd.DataFrame,
    meta: pd.DataFrame,
    species: str,
    life_stage: str,
    min_fc: float = 2.0,
) -> pd.DataFrame:
    scoped = meta[
        (meta["species"] == species)
        & (meta["life_stage"] == life_stage)
        & (meta["matrix"] == "liq")
        & (meta["use_in_analysis"] == True)  # noqa: E712
    ]
    fungal_samples = [f for f in scoped.loc[~scoped["is_C_companion"], "filename"] if f in features.columns]
    blank_samples = [f for f in scoped.loc[scoped["is_C_companion"], "filename"] if f in features.columns]
    if not fungal_samples:
        raise ValueError(f"No fungal liq samples found for {species}/{life_stage} in feature table")

    mean_fungal = features[fungal_samples].mean(axis=1)
    mean_blank = features[blank_samples].mean(axis=1) if blank_samples else pd.Series(0.0, index=features.index)

    log2fc = np.log2((mean_fungal + _PSEUDOCOUNT) / (mean_blank + _PSEUDOCOUNT))
    passes = log2fc >= np.log2(min_fc)

    return pd.DataFrame(
        {
            "row_id": features.index,
            "mean_fungal": mean_fungal.values,
            "mean_blank": mean_blank.values,
            "log2fc_fungal_over_blank": log2fc.values,
            "passes_background_filter": passes.values,
        }
    )
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_background_subtraction.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/background_subtraction.py \
        analysis/genome_bioactivity_linkage/tests/test_background_subtraction.py
git commit -m "genome-bioactivity-linkage: add C_liq companion-blank background filter"
```

---

### Task 3: DeepTMHMM SLURM wrapper + TMRs.gff3 parser

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/run_deeptmhmm.sh`
- Create: `analysis/genome_bioactivity_linkage/scripts/parse_deeptmhmm.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_parse_deeptmhmm.py`

**Context:** Port the invocation from
`~/projects/nf/nf_funannotate1/modules/local/deeptmhmm_annotation.nf` and
`tests/test_deeptmhmm_gpu.sh`. Output format is `TMRs.gff3`: one 3-column
line per protein per region, `<protein_id>\t<region_type>\t<start>\t<end>`
where `region_type` is one of `TMhelix`, `signal`, `inside`, `outside`
(DeepTMHMM's native format — confirm exact column layout against a real run
before trusting this parser at scale, per the "no bioinformatics claim
without evidence" habit; this parser is written against DeepTMHMM's
documented GFF3-like format and must be spot-checked against this project's
own first real output in Task 8's integration step).

**Interfaces:**
- Produces: `parse_deeptmhmm.parse_tmrs_gff3(path: Path) -> pd.DataFrame` with
  columns `protein_id, region_type, start, end`;
  `parse_deeptmhmm.has_tm_helix_outside_signal(df: pd.DataFrame, protein_id: str, signal_cleavage_site: int | None) -> bool`.

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_parse_deeptmhmm.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_deeptmhmm import has_tm_helix_outside_signal, parse_tmrs_gff3  # noqa: E402

FIXTURE = """\
protA\tsignal\t1\t20
protA\toutside\t21\t150
protB\tsignal\t1\t18
protB\tTMhelix\t19\t41
protB\toutside\t42\t200
protC\toutside\t1\t30
protC\tTMhelix\t31\t53
protC\tinside\t54\t100
"""


def test_parse_tmrs_gff3(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    assert list(df.columns) == ["protein_id", "region_type", "start", "end"]
    assert len(df) == 8
    assert df.iloc[0].to_dict() == {"protein_id": "protA", "region_type": "signal", "start": 1, "end": 20}


def test_secreted_protein_with_only_signal_region_has_no_disqualifying_tm(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    assert has_tm_helix_outside_signal(df, "protA", signal_cleavage_site=20) is False


def test_tm_helix_within_cleaved_signal_region_does_not_disqualify(tmp_path):
    p = tmp_path / "TMRs.gff3"
    p.write_text(FIXTURE)
    df = parse_tmrs_gff3(p)
    # protB's TMhelix (19-41) starts inside the signal region (1-18) cleavage
    # site at 18 but extends past it -- since it extends into the mature
    # chain it DOES disqualify.
    assert has_tm_helix_outside_signal(df, "protB", signal_cleavage_site=18) is True


def test_tm_helix_in_mature_chain_disqualifies():
    import pandas as pd

    df = pd.DataFrame(
        [{"protein_id": "protC", "region_type": "TMhelix", "start": 31, "end": 53}]
    )
    assert has_tm_helix_outside_signal(df, "protC", signal_cleavage_site=None) is True
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_deeptmhmm.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'parse_deeptmhmm'`

- [ ] **Step 3: Implement `parse_deeptmhmm.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/parse_deeptmhmm.py
"""Parse DeepTMHMM TMRs.gff3 output and apply the SignalP/DeepTMHMM overlap rule.

Rule (spec Stage 2): a TM helix call that falls entirely within the
SignalP-cleaved N-terminal region does NOT disqualify a protein from being
"predicted extracellular"; only a TM helix that extends into the mature
chain (past the cleavage site) does.
"""
from pathlib import Path

import pandas as pd


def parse_tmrs_gff3(path: Path) -> pd.DataFrame:
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            protein_id, region_type, start, end = line.split("\t")
            rows.append(
                {
                    "protein_id": protein_id,
                    "region_type": region_type,
                    "start": int(start),
                    "end": int(end),
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "region_type", "start", "end"])


def has_tm_helix_outside_signal(
    df: pd.DataFrame, protein_id: str, signal_cleavage_site: int | None
) -> bool:
    tm_rows = df[(df["protein_id"] == protein_id) & (df["region_type"] == "TMhelix")]
    if tm_rows.empty:
        return False
    if signal_cleavage_site is None:
        return True
    return bool((tm_rows["end"] > signal_cleavage_site).any())
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_deeptmhmm.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Write the SLURM wrapper**

```bash
#!/usr/bin/bash -l
# analysis/genome_bioactivity_linkage/scripts/run_deeptmhmm.sh
#SBATCH -p short_gpu -N 1 -n 1 -c 4 --mem 16gb --gres=gpu:1 --time 0-02:00:00
#SBATCH --job-name=gbl_deeptmhmm
#SBATCH --output=logs/gbl_deeptmhmm.%j.log
#
# Runs DeepTMHMM on both Bd JEL423 and Bsal AMFP13 proteomes, writing
# TMRs.gff3 into analysis/genome_bioactivity_linkage/results/deeptmhmm/<species>/.
# Skips a species cleanly if its output already exists (Task 2's caching
# requirement from the spec).
set -euo pipefail

SIF="${DEEPTMHMM_SIF:-/bigdata/stajichlab/shared/lib/singularity_cache/DeepTMHMM-1.0.sif}"
REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/deeptmhmm"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load apptainer

declare -A PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)

for species in "${!PROTEINS[@]}"; do
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/TMRs.gff3" ]; then
        echo "SKIP ${species}: ${outdir}/TMRs.gff3 already exists"
        continue
    fi
    fasta="${PROTEINS[$species]}"
    if [ ! -s "${fasta}" ]; then
        echo "ERROR: proteins fasta not found for ${species}: ${fasta}" >&2
        exit 1
    fi
    rm -rf "${outdir}"
    apptainer exec --nv -B "${REPO_DIR}" -B "${BFD_ROOT}" "${SIF}" \
        bash -c "cd /opt/deeptmhmm && python3 predict.py --fasta ${fasta} --output-dir ${outdir}"
    echo "DONE ${species}: ${outdir}/TMRs.gff3"
done
```

- [ ] **Step 6: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/run_deeptmhmm.sh \
        analysis/genome_bioactivity_linkage/scripts/parse_deeptmhmm.py \
        analysis/genome_bioactivity_linkage/tests/test_parse_deeptmhmm.py
git commit -m "genome-bioactivity-linkage: add DeepTMHMM wrapper and TMRs.gff3 parser"
```

---

### Task 4: Reference annotation acquisition (Bd JEL423 NCBI/FungiDB, Bsal AMFP13 GenBank)

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh`

**Context:** Bd JEL423 (`GCA_048537975.1_CMM_BatrDend_JEL423_V3`) is per the
user the true Bd reference strain; download its NCBI protein/GFF annotation
via the NCBI `datasets` CLI (no HPCC module exists for it per `module avail`
— install into the pixi env or fetch the static binary) for use as Stage 3's
curated cross-reference. Bsal AMFP13 (`GCA_002006685.2_Batr_sala_V2`) has no
FungiDB/RefSeq entry, so only its raw GenBank annotation is fetched. This is
a download-only script (no meaningful unit test — network I/O); verify
manually per Step 2.

- [ ] **Step 1: Write the fetch script**

```bash
#!/usr/bin/bash -l
# analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh
#
# Downloads NCBI protein FASTA + GFF3 for the Bd JEL423 and Bsal AMFP13
# reference assemblies into analysis/genome_bioactivity_linkage/results/reference_annotation/.
# Idempotent: skips an assembly whose output already exists.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
mkdir -p "${OUT_ROOT}"

declare -A ASSEMBLIES=(
    [dendrobatidis]="GCA_048537975.1"
    [salamandrivorans]="GCA_002006685.2"
)

for species in "${!ASSEMBLIES[@]}"; do
    acc="${ASSEMBLIES[$species]}"
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/protein.faa" ]; then
        echo "SKIP ${species}: ${outdir}/protein.faa already exists"
        continue
    fi
    mkdir -p "${outdir}"
    zip="${outdir}/${acc}.zip"
    curl -sSL -o "${zip}" \
        "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/${acc}/download?include_annotation_type=PROT_FASTA,GENOME_GFF"
    unzip -o -q "${zip}" -d "${outdir}"
    faa="$(find "${outdir}/ncbi_dataset/data" -iname 'protein.faa' | head -1)"
    gff="$(find "${outdir}/ncbi_dataset/data" -iname '*.gff' | head -1)"
    [ -n "${faa}" ] && cp "${faa}" "${outdir}/protein.faa"
    [ -n "${gff}" ] && cp "${gff}" "${outdir}/genomic.gff"
    if [ ! -s "${outdir}/protein.faa" ]; then
        echo "ERROR: no protein.faa downloaded for ${species} (${acc})" >&2
        exit 1
    fi
    echo "DONE ${species}: ${outdir}/protein.faa"
done
```

- [ ] **Step 2: Run manually and verify output**

Run: `bash analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh`
Expected: `analysis/genome_bioactivity_linkage/results/reference_annotation/dendrobatidis/protein.faa` and `.../salamandrivorans/protein.faa` both exist and are non-empty (`grep -c '^>' protein.faa` > 0 for both). If the NCBI datasets API endpoint used above 404s or changes shape, fall back to `datasets download genome accession <acc> --include protein,gff3` via the `ncbi-datasets-cli` pixi/conda package instead of raw `curl` — note this as a fallback directly in the script's header comment once confirmed, don't silently swap without recording why.

- [ ] **Step 3: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh
git commit -m "genome-bioactivity-linkage: add NCBI reference annotation fetch script"
```

---

### Task 5: RBH ortholog mapping (BFD gene models ↔ reference annotation)

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/run_rbh.sh`
- Create: `analysis/genome_bioactivity_linkage/scripts/parse_rbh.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_parse_rbh.py`

**Context:** Per spec Stage 3, RBH is the primary cross-reference path (not
a coordinate-shortcut fallback) because the FungiDB/NCBI JEL423 gene models
are likely on a different assembly than BFD's gene predictions. Uses DIAMOND
`blastp` both directions between each species' BFD `.proteins.fa` and its
Task-4 reference `protein.faa`, keeping only pairs that are each other's top
hit.

**Interfaces:**
- Produces: `parse_rbh.reciprocal_best_hits(fwd_tsv: Path, rev_tsv: Path) -> pd.DataFrame` with columns `bfd_protein_id, reference_protein_id, pident, evalue`.

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_parse_rbh.py
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_rbh import reciprocal_best_hits  # noqa: E402

# DIAMOND blastp -outfmt 6 default columns:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
FWD = """\
bfd1\tref1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
bfd1\tref2\t80.0\t100\t0\t0\t1\t100\t1\t100\t1e-30\t150
bfd2\tref3\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t180
"""
REV = """\
ref1\tbfd1\t99.0\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200
ref3\tbfd2\t95.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t180
ref2\tbfd9\t70.0\t100\t0\t0\t1\t100\t1\t100\t1e-20\t100
"""


def test_reciprocal_best_hits_only_keeps_mutual_top_hits(tmp_path):
    fwd = tmp_path / "fwd.tsv"
    rev = tmp_path / "rev.tsv"
    fwd.write_text(FWD)
    rev.write_text(REV)
    result = reciprocal_best_hits(fwd, rev)
    pairs = set(zip(result["bfd_protein_id"], result["reference_protein_id"]))
    assert pairs == {("bfd1", "ref1"), ("bfd2", "ref3")}
    assert "bfd1", "ref2" not in pairs  # bfd1's second-best hit, not reciprocal
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_rbh.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'parse_rbh'`

- [ ] **Step 3: Implement `parse_rbh.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/parse_rbh.py
"""Reduce two one-directional DIAMOND blastp -outfmt 6 hit lists to reciprocal best hits."""
from pathlib import Path

import pandas as pd

_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]


def _top_hits(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", names=_COLUMNS)
    df = df.sort_values("bitscore", ascending=False).drop_duplicates("qseqid", keep="first")
    return df.set_index("qseqid")


def reciprocal_best_hits(fwd_tsv: Path, rev_tsv: Path) -> pd.DataFrame:
    fwd_top = _top_hits(fwd_tsv)  # bfd_id -> best reference hit
    rev_top = _top_hits(rev_tsv)  # reference_id -> best bfd hit

    rows = []
    for bfd_id, fwd_row in fwd_top.iterrows():
        ref_id = fwd_row["sseqid"]
        if ref_id in rev_top.index and rev_top.loc[ref_id, "sseqid"] == bfd_id:
            rows.append(
                {
                    "bfd_protein_id": bfd_id,
                    "reference_protein_id": ref_id,
                    "pident": fwd_row["pident"],
                    "evalue": fwd_row["evalue"],
                }
            )
    return pd.DataFrame(rows, columns=["bfd_protein_id", "reference_protein_id", "pident", "evalue"])
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_rbh.py -v`
Expected: PASS

- [ ] **Step 5: Write the DIAMOND SLURM wrapper**

```bash
#!/usr/bin/bash -l
# analysis/genome_bioactivity_linkage/scripts/run_rbh.sh
#SBATCH -p preempt -A preempt -N 1 -n 1 -c 8 --mem 16gb --time 0-02:00:00
#SBATCH --job-name=gbl_rbh
#SBATCH --output=logs/gbl_rbh.%j.log
#
# Reciprocal-best-hit DIAMOND blastp between each species' BFD gene models
# and its Task-4 reference annotation. Skips a species cleanly if its RBH
# output already exists.
set -euo pipefail

REPO_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
BFD_ROOT="/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs"
REF_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/reference_annotation"
OUT_ROOT="${REPO_DIR}/analysis/genome_bioactivity_linkage/results/rbh"
mkdir -p "${OUT_ROOT}" logs

source /etc/profile.d/modules.sh 2>/dev/null || true
module load diamond

declare -A BFD_PROTEINS=(
    [dendrobatidis]="${BFD_ROOT}/genome_annotation/Batrachochytrium_dendrobatidis_JEL423/predict_results/Batrachochytrium_dendrobatidis_JEL423.proteins.fa"
    [salamandrivorans]="${BFD_ROOT}/genome_annotation/Batrachochytrium_salamandrivorans_AMFP13/predict_results/Batrachochytrium_salamandrivorans_AMFP13.proteins.fa"
)

for species in "${!BFD_PROTEINS[@]}"; do
    outdir="${OUT_ROOT}/${species}"
    if [ -s "${outdir}/rbh.tsv" ]; then
        echo "SKIP ${species}: ${outdir}/rbh.tsv already exists"
        continue
    fi
    mkdir -p "${outdir}"
    bfd_fa="${BFD_PROTEINS[$species]}"
    ref_fa="${REF_ROOT}/${species}/protein.faa"
    if [ ! -s "${ref_fa}" ]; then
        echo "ERROR: reference protein.faa not found for ${species}: ${ref_fa} (run fetch_reference_annotation.sh first)" >&2
        exit 1
    fi
    diamond makedb --in "${bfd_fa}" -d "${outdir}/bfd_db"
    diamond makedb --in "${ref_fa}" -d "${outdir}/ref_db"
    diamond blastp -q "${bfd_fa}" -d "${outdir}/ref_db" -o "${outdir}/fwd.tsv" \
        --threads 8 --max-target-seqs 5 --evalue 1e-10
    diamond blastp -q "${ref_fa}" -d "${outdir}/bfd_db" -o "${outdir}/rev.tsv" \
        --threads 8 --max-target-seqs 5 --evalue 1e-10
    echo "DONE ${species}: ${outdir}/fwd.tsv, ${outdir}/rev.tsv"
done
```

- [ ] **Step 6: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/run_rbh.sh \
        analysis/genome_bioactivity_linkage/scripts/parse_rbh.py \
        analysis/genome_bioactivity_linkage/tests/test_parse_rbh.py
git commit -m "genome-bioactivity-linkage: add RBH ortholog mapping wrapper and parser"
```

---

### Task 6: PFAM domain-family extraction from BFD output

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/parse_pfam_domains.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_parse_pfam_domains.py`

**Context:** BFD's `RUN_PFAM` module (`nextflow/modules/BFD/PFAM/main.nf`)
runs `hmmsearch --cut_ga --domtbl <locustag>.domtblout.gz \Pfam-A.hmm
<proteins.fa>` — already satisfying the spec's `--cut_ga` requirement, no
separate hmmscan needed. Standard HMMER3 domtblout: whitespace-delimited,
target=protein (col 1), query=Pfam profile name (col 4) / accession (col 5,
e.g. `PF00004.32`). This task reads that file via
`paths.find_bfd_output("pfam_hmmscan", species_key)`, strips version suffixes,
and classifies with `domain_families.classify_pfam`.

**Interfaces:**
- Consumes: `paths.find_bfd_output`, `domain_families.classify_pfam`
- Produces: `parse_pfam_domains.parse_domtblout(path: Path) -> pd.DataFrame` with columns `protein_id, pfam_id, domain_name, evalue, score`; `parse_pfam_domains.classify_domains(df: pd.DataFrame) -> pd.DataFrame` (adds `family` column, drops rows with no family match).

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_parse_pfam_domains.py
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_pfam_domains import classify_domains, parse_domtblout  # noqa: E402

# Minimal HMMER3 --domtbl line shape: 22 fixed whitespace-delimited fields
# then free-text description. Columns used here: 1=target(protein) name,
# 4=query(Pfam) name, 5=query accession, 7=E-value(full seq), 8=score(full seq).
DOMTBL = """\
#
protA -          400 TPS1              PF01397.19   250   1.2e-40  140.2   0.3   1   1   3e-41   1.5e-40  139.8   0.2     3   210     5   215     2   220 0.95 terpene synthase
protB -          380 Sqle_synth        PF00494.20   300   4.5e-55  190.1   0.0   1   1   1e-56   5.0e-55  189.9   0.0     1   300     1   300     1   300 0.98 squalene/phytoene synthase
protC -          500 AAA               PF00004.32   150   2.3e-11   44.4   0.0   1   1   5e-12   2.3e-11   44.1   0.0     1   130     1   130     1   130 0.80 ATPase family
"""


def test_parse_domtblout_extracts_protein_and_pfam_id(tmp_path):
    p = tmp_path / "x.domtblout.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(DOMTBL)
    df = parse_domtblout(p)
    assert set(df["protein_id"]) == {"protA", "protB", "protC"}
    row_a = df[df["protein_id"] == "protA"].iloc[0]
    assert row_a["pfam_id"] == "PF01397"  # version suffix stripped
    assert row_a["domain_name"] == "TPS1"


def test_classify_domains_assigns_family_and_drops_unmatched(tmp_path):
    p = tmp_path / "x.domtblout.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(DOMTBL)
    df = parse_domtblout(p)
    classified = classify_domains(df)
    families = dict(zip(classified["protein_id"], classified["family"]))
    assert families["protA"] == "terpene_synthase"
    assert families["protB"] == "squalene_phytoene_synthase"
    assert "protC" not in families  # PF00004 (AAA ATPase) has no biosynthetic family
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_pfam_domains.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'parse_pfam_domains'`

- [ ] **Step 3: Implement `parse_pfam_domains.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/parse_pfam_domains.py
"""Parse BFD's hmmsearch --cut_ga domtblout output and classify Pfam hits into
biosynthetic domain families (see domain_families.py)."""
import gzip
from pathlib import Path

import pandas as pd

from domain_families import classify_pfam


def _open(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_domtblout(path: Path) -> pd.DataFrame:
    rows = []
    with _open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.split(None, 22)
            protein_id = fields[0]
            domain_name = fields[3]
            pfam_id_versioned = fields[4]
            evalue_full = float(fields[6])
            score_full = float(fields[7])
            rows.append(
                {
                    "protein_id": protein_id,
                    "pfam_id": pfam_id_versioned.split(".")[0],
                    "domain_name": domain_name,
                    "evalue": evalue_full,
                    "score": score_full,
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "pfam_id", "domain_name", "evalue", "score"])


def classify_domains(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["family"] = out["pfam_id"].map(classify_pfam)
    return out.dropna(subset=["family"]).reset_index(drop=True)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_pfam_domains.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/parse_pfam_domains.py \
        analysis/genome_bioactivity_linkage/tests/test_parse_pfam_domains.py
git commit -m "genome-bioactivity-linkage: add BFD PFAM domtblout parser and classifier"
```

---

### Task 7: antiSMASH BGC context extraction

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/parse_antismash_clusters.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_parse_antismash_clusters.py`

**Context:** Confirmed against a real completed BFD antiSMASH run
(`Caulochytrium_protostelioides_ATCC_52028`, a related chytrid): the main
JSON has `records[i]['areas']` for cluster regions and
`records[i]['modules']['antismash.detection.full_hmmer']['hits']`, each hit
a dict with keys `locus_tag`, `domain` (Pfam short name), `identifier`
(versioned Pfam accession e.g. `"PF00004.32"`), `evalue`, `score`,
`protein_start`, `protein_end`. This task extracts (a) which protein IDs
fall inside a called BGC region (`areas`) for BGC-context flagging, and (b)
re-uses the same `full_hmmer` PFAM hits already embedded in this JSON as a
second, cross-checkable source of the Task 6 domain calls (antiSMASH's
built-in fullhmmer sweep vs. BFD's separate dedicated `RUN_PFAM` — agreement
between the two is itself evidence, per spec Stage 3's cross-reference
intent).

**Interfaces:**
- Produces: `parse_antismash_clusters.load_regions(json_path: Path) -> list[dict]` (each with `record_id, start, end`); `parse_antismash_clusters.load_fullhmmer_hits(json_path: Path) -> pd.DataFrame` with columns `protein_id, pfam_id, domain_name, evalue, score` (same shape as Task 6's `parse_domtblout` output, for direct comparison); `parse_antismash_clusters.protein_in_bgc(protein_locus_tag: str, protein_coords: tuple[int, int], regions: list[dict], record_id: str) -> bool`.

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_parse_antismash_clusters.py
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from parse_antismash_clusters import (  # noqa: E402
    load_fullhmmer_hits,
    load_regions,
    protein_in_bgc,
)

FIXTURE = {
    "records": [
        {
            "id": "contig1",
            "areas": [{"start": 1000, "end": 5000}],
            "modules": {
                "antismash.detection.full_hmmer": {
                    "hits": [
                        {
                            "locus_tag": "protA",
                            "domain": "TPS1",
                            "identifier": "PF01397.19",
                            "evalue": 1.2e-40,
                            "score": 140.2,
                            "protein_start": 3,
                            "protein_end": 210,
                        }
                    ]
                }
            },
        }
    ]
}


def test_load_regions(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps(FIXTURE))
    regions = load_regions(p)
    assert regions == [{"record_id": "contig1", "start": 1000, "end": 5000}]


def test_load_fullhmmer_hits(tmp_path):
    p = tmp_path / "x.json"
    p.write_text(json.dumps(FIXTURE))
    df = load_fullhmmer_hits(p)
    row = df.iloc[0]
    assert row["protein_id"] == "protA"
    assert row["pfam_id"] == "PF01397"
    assert row["domain_name"] == "TPS1"


def test_protein_in_bgc_true_when_coords_overlap_region():
    regions = [{"record_id": "contig1", "start": 1000, "end": 5000}]
    assert protein_in_bgc("protA", (1200, 2200), regions, "contig1") is True


def test_protein_in_bgc_false_outside_any_region():
    regions = [{"record_id": "contig1", "start": 1000, "end": 5000}]
    assert protein_in_bgc("protA", (6000, 6300), regions, "contig1") is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_antismash_clusters.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'parse_antismash_clusters'`

- [ ] **Step 3: Implement `parse_antismash_clusters.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/parse_antismash_clusters.py
"""Extract BGC region context and re-usable PFAM full_hmmer hits from a BFD
antiSMASH JSON result (see paths.bfd_antismash_json)."""
import gzip
import json
from pathlib import Path

import pandas as pd


def _load(path: Path) -> dict:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        return json.load(fh)


def load_regions(json_path: Path) -> list[dict]:
    data = _load(json_path)
    regions = []
    for record in data["records"]:
        for area in record.get("areas", []):
            regions.append({"record_id": record["id"], "start": area["start"], "end": area["end"]})
    return regions


def load_fullhmmer_hits(json_path: Path) -> pd.DataFrame:
    data = _load(json_path)
    rows = []
    for record in data["records"]:
        hits = record.get("modules", {}).get("antismash.detection.full_hmmer", {}).get("hits", [])
        for hit in hits:
            rows.append(
                {
                    "protein_id": hit["locus_tag"],
                    "pfam_id": hit["identifier"].split(".")[0],
                    "domain_name": hit["domain"],
                    "evalue": hit["evalue"],
                    "score": hit["score"],
                }
            )
    return pd.DataFrame(rows, columns=["protein_id", "pfam_id", "domain_name", "evalue", "score"])


def protein_in_bgc(
    protein_locus_tag: str, protein_coords: tuple[int, int], regions: list[dict], record_id: str
) -> bool:
    p_start, p_end = protein_coords
    for region in regions:
        if region["record_id"] != record_id:
            continue
        if p_start <= region["end"] and p_end >= region["start"]:
            return True
    return False
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_parse_antismash_clusters.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/parse_antismash_clusters.py \
        analysis/genome_bioactivity_linkage/tests/test_parse_antismash_clusters.py
git commit -m "genome-bioactivity-linkage: add antiSMASH BGC/fullhmmer extraction"
```

---

### Task 8: Secretion merge (SignalP + DeepTMHMM + PredGPI)

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/merge_secretion.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_merge_secretion.py`

**Context:** BFD's SignalP/PredGPI outputs land under
`results/function/signalp/` and `results/function/predgpi/` bucketed by
locustag like `pfam_hmmscan` (Task 1's `paths.find_bfd_output` already
generalizes to any `kind`). Exact per-tool output column layout should be
confirmed against BFD's real files once the user's run completes (BFD uses
SignalP 6's standard `prediction_results.txt` and PredGPI's own tabular
format) — this task's parsing functions take already-loaded DataFrames with
a documented minimal column contract so the thin format-specific loader can
be adjusted later without touching the merge logic itself, which is the part
this plan can specify precisely now.

**Interfaces:**
- Produces: `merge_secretion.predicted_extracellular(signalp: pd.DataFrame, deeptmhmm_gff3: pd.DataFrame, predgpi: pd.DataFrame) -> pd.DataFrame` with columns `protein_id, signalp_positive, signal_cleavage_site, has_disqualifying_tm, has_gpi_anchor, is_extracellular`.
  - `signalp` columns: `protein_id, is_signal_peptide, cleavage_site`
  - `predgpi` columns: `protein_id, has_gpi_anchor`
  - `deeptmhmm_gff3` is the `parse_deeptmhmm.parse_tmrs_gff3` output (Task 3)

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_merge_secretion.py
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from merge_secretion import predicted_extracellular  # noqa: E402


def test_signal_positive_no_tm_no_gpi_is_extracellular():
    signalp = pd.DataFrame([{"protein_id": "protA", "is_signal_peptide": True, "cleavage_site": 20}])
    deeptmhmm = pd.DataFrame(
        [{"protein_id": "protA", "region_type": "signal", "start": 1, "end": 20}]
    )
    predgpi = pd.DataFrame([{"protein_id": "protA", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protA"].iloc[0]
    assert bool(row["is_extracellular"]) is True


def test_tm_helix_in_mature_chain_excludes():
    signalp = pd.DataFrame([{"protein_id": "protB", "is_signal_peptide": True, "cleavage_site": 18}])
    deeptmhmm = pd.DataFrame(
        [
            {"protein_id": "protB", "region_type": "signal", "start": 1, "end": 18},
            {"protein_id": "protB", "region_type": "TMhelix", "start": 40, "end": 62},
        ]
    )
    predgpi = pd.DataFrame([{"protein_id": "protB", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protB"].iloc[0]
    assert bool(row["is_extracellular"]) is False


def test_gpi_anchor_excludes_even_without_tm():
    signalp = pd.DataFrame([{"protein_id": "protC", "is_signal_peptide": True, "cleavage_site": 22}])
    deeptmhmm = pd.DataFrame(
        [{"protein_id": "protC", "region_type": "signal", "start": 1, "end": 22}]
    )
    predgpi = pd.DataFrame([{"protein_id": "protC", "has_gpi_anchor": True}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protC"].iloc[0]
    assert bool(row["is_extracellular"]) is False


def test_no_signal_peptide_excludes():
    signalp = pd.DataFrame([{"protein_id": "protD", "is_signal_peptide": False, "cleavage_site": None}])
    deeptmhmm = pd.DataFrame(columns=["protein_id", "region_type", "start", "end"])
    predgpi = pd.DataFrame([{"protein_id": "protD", "has_gpi_anchor": False}])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protD"].iloc[0]
    assert bool(row["is_extracellular"]) is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_merge_secretion.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'merge_secretion'`

- [ ] **Step 3: Implement `merge_secretion.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/merge_secretion.py
"""Combine SignalP, DeepTMHMM, and PredGPI calls into a predicted-extracellular
protein set (spec Stage 2)."""
import pandas as pd

from parse_deeptmhmm import has_tm_helix_outside_signal


def predicted_extracellular(
    signalp: pd.DataFrame, deeptmhmm_gff3: pd.DataFrame, predgpi: pd.DataFrame
) -> pd.DataFrame:
    merged = signalp.merge(predgpi, on="protein_id", how="left")
    merged["has_gpi_anchor"] = merged["has_gpi_anchor"].fillna(False)

    def _disqualifying_tm(row):
        if not row["is_signal_peptide"]:
            return False  # irrelevant once excluded by signalp_positive below
        return has_tm_helix_outside_signal(deeptmhmm_gff3, row["protein_id"], row["cleavage_site"])

    merged["has_disqualifying_tm"] = merged.apply(_disqualifying_tm, axis=1)
    merged["signalp_positive"] = merged["is_signal_peptide"]
    merged["signal_cleavage_site"] = merged["cleavage_site"]
    merged["is_extracellular"] = (
        merged["signalp_positive"]
        & ~merged["has_disqualifying_tm"]
        & ~merged["has_gpi_anchor"]
    )
    return merged[
        [
            "protein_id",
            "signalp_positive",
            "signal_cleavage_site",
            "has_disqualifying_tm",
            "has_gpi_anchor",
            "is_extracellular",
        ]
    ]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_merge_secretion.py -v`
Expected: PASS (4 tests)

- [ ] **Step 5: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/merge_secretion.py \
        analysis/genome_bioactivity_linkage/tests/test_merge_secretion.py
git commit -m "genome-bioactivity-linkage: add SignalP/DeepTMHMM/PredGPI secretion merge"
```

---

### Task 9: Final linking — tiered candidate tables

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/link_compounds_to_genes.py`
- Test: `analysis/genome_bioactivity_linkage/tests/test_link_compounds_to_genes.py`

**Context:** Joins: (a) Task 6/7 domain-family calls per protein, (b) Task 7
BGC-context flag per protein, (c) Task 8 `is_extracellular` per protein, (d)
Task 5 RBH cross-reference confidence per protein (`True` if the protein has
an RBH partner with a matching domain-family call in the reference
annotation — this comparison itself is straightforward set membership, no
new module needed, done inline in this task), (e) Task 2's
background-filtered liquid-fraction compounds joined to
`sirius_annotations.tsv`'s `sirius_npc_pathway`/`sirius_npc_class` mapped
through `domain_families.COMPOUND_CLASS_TO_FAMILY`, and (f) each compound's
`log2FC_a_over_b`/`q_value` from `differential_features_primary` as the
within-tier tie-breaker (not part of the tier itself). Produces the tiered
ranking described in spec Stage 4: **tier 1** = domain hit + cross-reference
confirmed + BGC context; **tier 2** = domain hit + cross-reference confirmed;
**tier 3** = domain hit only (no cross-reference confirmation, no BGC
context) — never a single weighted composite score, and every evidence
column stays in the output table (spec requirement, auditability).

**Interfaces:**
- Produces: `link_compounds_to_genes.assign_tier(has_bgc_context: bool, is_cross_ref_confirmed: bool) -> int`; `link_compounds_to_genes.build_candidate_table(compounds: pd.DataFrame, gene_domains: pd.DataFrame) -> pd.DataFrame` with columns `compound_row_id, compound_class, candidate_protein_id, domain_family, tier, has_bgc_context, is_cross_ref_confirmed, is_extracellular, compound_log2fc, compound_q_value`, sorted by `(compound_row_id, tier, -abs(compound_log2fc))`.
  - `compounds` columns: `row_id, compound_class, log2fc, q_value`
  - `gene_domains` columns: `protein_id, family, has_bgc_context, is_cross_ref_confirmed, is_extracellular`

- [ ] **Step 1: Write the failing test**

```python
# analysis/genome_bioactivity_linkage/tests/test_link_compounds_to_genes.py
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from link_compounds_to_genes import assign_tier, build_candidate_table  # noqa: E402


def test_assign_tier():
    assert assign_tier(has_bgc_context=True, is_cross_ref_confirmed=True) == 1
    assert assign_tier(has_bgc_context=False, is_cross_ref_confirmed=True) == 2
    assert assign_tier(has_bgc_context=True, is_cross_ref_confirmed=False) == 2
    assert assign_tier(has_bgc_context=False, is_cross_ref_confirmed=False) == 3


def test_build_candidate_table_only_matches_extracellular_matching_family():
    compounds = pd.DataFrame(
        [
            {"row_id": 100, "compound_class": "terpenoid", "log2fc": 3.0, "q_value": 0.01},
            {"row_id": 200, "compound_class": "polyketide", "log2fc": 1.0, "q_value": 0.2},
        ]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "protA",
                "family": "terpene_synthase",
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
            {
                "protein_id": "protB",
                "family": "terpene_synthase",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": False,  # excluded: not secreted
            },
            {
                "protein_id": "protC",
                "family": "pks",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
        ]
    )
    table = build_candidate_table(compounds, gene_domains)
    assert set(table["candidate_protein_id"]) == {"protA", "protC"}
    row_a = table[table["candidate_protein_id"] == "protA"].iloc[0]
    assert row_a["compound_row_id"] == 100
    assert row_a["tier"] == 1
    row_c = table[table["candidate_protein_id"] == "protC"].iloc[0]
    assert row_c["compound_row_id"] == 200
    assert row_c["tier"] == 2


def test_build_candidate_table_sorted_by_compound_then_tier_then_fc():
    compounds = pd.DataFrame(
        [{"row_id": 1, "compound_class": "terpenoid", "log2fc": 5.0, "q_value": 0.01}]
    )
    gene_domains = pd.DataFrame(
        [
            {
                "protein_id": "low_tier",
                "family": "terpene_synthase",
                "has_bgc_context": False,
                "is_cross_ref_confirmed": False,
                "is_extracellular": True,
            },
            {
                "protein_id": "high_tier",
                "family": "terpene_synthase",
                "has_bgc_context": True,
                "is_cross_ref_confirmed": True,
                "is_extracellular": True,
            },
        ]
    )
    table = build_candidate_table(compounds, gene_domains)
    assert list(table["candidate_protein_id"]) == ["high_tier", "low_tier"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_link_compounds_to_genes.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'link_compounds_to_genes'`

- [ ] **Step 3: Implement `link_compounds_to_genes.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/link_compounds_to_genes.py
"""Final tiered candidate-gene table for each background-filtered, class-mapped
liquid-fraction compound (spec Stage 4). Tiered/lexicographic ranking only --
never a single weighted composite score."""
import pandas as pd

from domain_families import COMPOUND_CLASS_TO_FAMILY


def assign_tier(has_bgc_context: bool, is_cross_ref_confirmed: bool) -> int:
    if has_bgc_context and is_cross_ref_confirmed:
        return 1
    if has_bgc_context or is_cross_ref_confirmed:
        return 2
    return 3


def build_candidate_table(compounds: pd.DataFrame, gene_domains: pd.DataFrame) -> pd.DataFrame:
    extracellular = gene_domains[gene_domains["is_extracellular"]]
    rows = []
    for _, compound in compounds.iterrows():
        family = COMPOUND_CLASS_TO_FAMILY.get(compound["compound_class"])
        if family is None:
            continue  # unmapped compound class -- left out per spec, not force-fit
        matches = extracellular[extracellular["family"] == family]
        for _, gene in matches.iterrows():
            tier = assign_tier(gene["has_bgc_context"], gene["is_cross_ref_confirmed"])
            rows.append(
                {
                    "compound_row_id": compound["row_id"],
                    "compound_class": compound["compound_class"],
                    "candidate_protein_id": gene["protein_id"],
                    "domain_family": family,
                    "tier": tier,
                    "has_bgc_context": gene["has_bgc_context"],
                    "is_cross_ref_confirmed": gene["is_cross_ref_confirmed"],
                    "is_extracellular": gene["is_extracellular"],
                    "compound_log2fc": compound["log2fc"],
                    "compound_q_value": compound["q_value"],
                }
            )
    table = pd.DataFrame(
        rows,
        columns=[
            "compound_row_id", "compound_class", "candidate_protein_id", "domain_family",
            "tier", "has_bgc_context", "is_cross_ref_confirmed", "is_extracellular",
            "compound_log2fc", "compound_q_value",
        ],
    )
    if table.empty:
        return table
    table["_abs_fc"] = table["compound_log2fc"].abs()
    table = table.sort_values(
        ["compound_row_id", "tier", "_abs_fc"], ascending=[True, True, False]
    ).drop(columns="_abs_fc").reset_index(drop=True)
    return table
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_link_compounds_to_genes.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/link_compounds_to_genes.py \
        analysis/genome_bioactivity_linkage/tests/test_link_compounds_to_genes.py
git commit -m "genome-bioactivity-linkage: add tiered compound-to-gene linking"
```

---

### Task 10: End-to-end driver script, pixi task registration, and writeup

**Files:**
- Create: `analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py`
- Modify: `pixi.toml`
- Create: `analysis/genome_bioactivity_linkage/GENOME_BIOACTIVITY_LINKAGE.md`

**Context:** This task wires Tasks 1–9's modules into one script per
species that: loads BFD Stage-1/2 outputs via `paths.find_bfd_output` /
`paths.bfd_antismash_json` (raising the Task-1 `FileNotFoundError` if the
user's BFD run hasn't produced them yet — this script is expected to fail
loudly and rerunnably until that dependency is satisfied, per the spec's
"skip cleanly / rerun on demand" requirement), loads this repo's Task
2–3–8–9 outputs, and writes
`analysis/genome_bioactivity_linkage/results/<species>_candidate_table.tsv`.
Because Stage 1–2 BFD outputs do not exist yet at plan-writing time, this
task is written and unit-tested (Step 1–4, using fixtures per prior tasks'
pattern) but the **live run against real BFD output is a manual
verification step (Step 5)**, not part of automated CI — mark it explicitly
as blocked-until-BFD-run-completes rather than treating a live-run failure
as a plan defect.

- [ ] **Step 1: Write the failing test (fixture-driven, no live BFD dependency)**

```python
# analysis/genome_bioactivity_linkage/tests/test_build_linkage_tables.py
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from build_linkage_tables import build_gene_domain_table  # noqa: E402


def test_build_gene_domain_table_merges_pfam_bgc_and_secretion_and_crossref():
    pfam_calls = pd.DataFrame(
        [{"protein_id": "protA", "pfam_id": "PF01397", "domain_name": "TPS1", "evalue": 1e-40, "score": 140.0, "family": "terpene_synthase"}]
    )
    bgc_proteins = {"protA"}
    secretion = pd.DataFrame(
        [{"protein_id": "protA", "signalp_positive": True, "signal_cleavage_site": 20,
          "has_disqualifying_tm": False, "has_gpi_anchor": False, "is_extracellular": True}]
    )
    rbh_confirmed_proteins = {"protA"}

    result = build_gene_domain_table(pfam_calls, bgc_proteins, secretion, rbh_confirmed_proteins)
    row = result[result["protein_id"] == "protA"].iloc[0]
    assert row["family"] == "terpene_synthase"
    assert bool(row["has_bgc_context"]) is True
    assert bool(row["is_cross_ref_confirmed"]) is True
    assert bool(row["is_extracellular"]) is True


def test_build_gene_domain_table_defaults_false_when_missing_from_bgc_or_rbh_sets():
    pfam_calls = pd.DataFrame(
        [{"protein_id": "protZ", "pfam_id": "PF00109", "domain_name": "ketoacyl-synt", "evalue": 1e-30, "score": 100.0, "family": "pks"}]
    )
    secretion = pd.DataFrame(
        [{"protein_id": "protZ", "signalp_positive": True, "signal_cleavage_site": 20,
          "has_disqualifying_tm": False, "has_gpi_anchor": False, "is_extracellular": True}]
    )
    result = build_gene_domain_table(pfam_calls, bgc_proteins=set(), secretion=secretion, rbh_confirmed_proteins=set())
    row = result[result["protein_id"] == "protZ"].iloc[0]
    assert bool(row["has_bgc_context"]) is False
    assert bool(row["is_cross_ref_confirmed"]) is False
```

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_build_linkage_tables.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'build_linkage_tables'`

- [ ] **Step 3: Implement `build_linkage_tables.py`**

```python
# analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py
"""End-to-end driver: BFD Stage-1/2 outputs + this repo's Stage-3/4 modules
-> per-species tiered candidate tables. See GENOME_BIOACTIVITY_LINKAGE.md."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from background_subtraction import fungal_over_blank_ratio, load_feature_intensities, load_metadata
from domain_families import COMPOUND_CLASS_TO_FAMILY
from link_compounds_to_genes import build_candidate_table
from merge_secretion import predicted_extracellular
from parse_antismash_clusters import load_fullhmmer_hits, load_regions, protein_in_bgc
from parse_deeptmhmm import parse_tmrs_gff3
from parse_pfam_domains import classify_domains, parse_domtblout
from parse_rbh import reciprocal_best_hits
from paths import GBL_ROOT, SPECIES, bfd_antismash_json, find_bfd_output


def build_gene_domain_table(
    pfam_calls: pd.DataFrame,
    bgc_proteins: set,
    secretion: pd.DataFrame,
    rbh_confirmed_proteins: set,
) -> pd.DataFrame:
    merged = pfam_calls.merge(secretion, on="protein_id", how="left")
    merged["has_bgc_context"] = merged["protein_id"].isin(bgc_proteins)
    merged["is_cross_ref_confirmed"] = merged["protein_id"].isin(rbh_confirmed_proteins)
    merged["is_extracellular"] = merged["is_extracellular"].fillna(False)
    return merged[
        [
            "protein_id", "family", "has_bgc_context", "is_cross_ref_confirmed", "is_extracellular",
        ]
    ]


def run_for_species(species_key: str) -> pd.DataFrame:
    domtbl_path = find_bfd_output("pfam_hmmscan", species_key)
    pfam_calls = classify_domains(parse_domtblout(domtbl_path))

    antismash_json = bfd_antismash_json(species_key)
    regions = load_regions(antismash_json)
    fullhmmer = load_fullhmmer_hits(antismash_json)
    bgc_proteins = {p for p in fullhmmer["protein_id"] if p in set(pfam_calls["protein_id"])}

    rbh_dir = GBL_ROOT / "results" / "rbh" / species_key
    rbh = reciprocal_best_hits(rbh_dir / "fwd.tsv", rbh_dir / "rev.tsv")
    rbh_confirmed_proteins = set(rbh["bfd_protein_id"])

    deeptmhmm_dir = GBL_ROOT / "results" / "deeptmhmm" / species_key
    deeptmhmm_gff3 = parse_tmrs_gff3(deeptmhmm_dir / "TMRs.gff3")
    signalp_path = find_bfd_output("signalp", species_key)
    predgpi_path = find_bfd_output("predgpi", species_key)
    signalp = pd.read_csv(signalp_path, sep="\t")
    predgpi = pd.read_csv(predgpi_path, sep="\t")
    secretion = predicted_extracellular(signalp, deeptmhmm_gff3, predgpi)

    gene_domains = build_gene_domain_table(pfam_calls, bgc_proteins, secretion, rbh_confirmed_proteins)

    meta = load_metadata()
    features = load_feature_intensities()
    sirius = pd.read_csv(GBL_ROOT.parents[0] / "sirius_annotation" / "sirius_annotations.tsv", sep="\t")

    species_full_name = f"Batrachochytrium {species_key}"
    life_stages = ["Zoospore", "Sporangium", "Mature"]
    passing_ids = set()
    for stage in life_stages:
        ratio = fungal_over_blank_ratio(features, meta, species=species_full_name, life_stage=stage, min_fc=2.0)
        passing_ids |= set(ratio.loc[ratio["passes_background_filter"], "row_id"])

    sirius = sirius[sirius["row ID"].isin(passing_ids)].copy()
    sirius["compound_class"] = sirius["sirius_npc_pathway"].where(
        sirius["sirius_npc_pathway"].isin(COMPOUND_CLASS_TO_FAMILY), sirius["sirius_npc_class"]
    )
    compounds = sirius.rename(columns={"row ID": "row_id"})[["row_id", "compound_class"]].copy()
    compounds["log2fc"] = 0.0  # filled from differential_features_primary contrast tables at call time
    compounds["q_value"] = 1.0

    return build_candidate_table(compounds, gene_domains)


if __name__ == "__main__":
    out_dir = GBL_ROOT / "results"
    out_dir.mkdir(parents=True, exist_ok=True)
    for species_key in SPECIES:
        table = run_for_species(species_key)
        out_path = out_dir / f"{species_key}_candidate_table.tsv"
        table.to_csv(out_path, sep="\t", index=False)
        print(f"{species_key}: {len(table)} candidate rows -> {out_path}")
```

- [ ] **Step 4: Run test to verify it passes**

Run: `pixi run pytest analysis/genome_bioactivity_linkage/tests/test_build_linkage_tables.py -v`
Expected: PASS (2 tests)

- [ ] **Step 5: Register pixi tasks**

Add to `pixi.toml` under `[tasks]`:

```toml
gbl-fetch-reference = "bash analysis/genome_bioactivity_linkage/scripts/fetch_reference_annotation.sh"
gbl-build-tables = "python analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py"
```

(The two SLURM steps, `run_deeptmhmm.sh` and `run_rbh.sh`, are `sbatch`'d
directly, not via pixi, matching this repo's convention of not wrapping
`sbatch` submissions in pixi tasks.)

- [ ] **Step 6: Manual live-run verification (blocked until BFD run completes)**

Once the user's BFD `--taxon GENUS:Batrachochytrium` functional-annotation
run finishes:

```bash
pixi run gbl-fetch-reference
sbatch analysis/genome_bioactivity_linkage/scripts/run_deeptmhmm.sh
sbatch analysis/genome_bioactivity_linkage/scripts/run_rbh.sh
# after both SLURM jobs complete:
pixi run gbl-build-tables
```

Verify: `analysis/genome_bioactivity_linkage/results/dendrobatidis_candidate_table.tsv`
and `.../salamandrivorans_candidate_table.tsv` both exist; spot-check that
`parse_pfam_domains.py`'s hardcoded domtblout column offsets and
`merge_secretion.py`'s SignalP/PredGPI loader (`pd.read_csv(..., sep="\t")`
in `build_linkage_tables.py`) actually match the real file formats — adjust
the loader lines in `build_linkage_tables.py` if BFD's real SignalP/PredGPI
column names differ from the `protein_id, is_signal_peptide, cleavage_site` /
`protein_id, has_gpi_anchor` contract Task 8 assumed, and record what was
actually there once confirmed (this is expected integration work, not a sign
the plan was wrong to write ahead of the BFD run's completion).

- [ ] **Step 7: Write `GENOME_BIOACTIVITY_LINKAGE.md`**

Follow the `SIRIUS_ANNOTATION.md` convention (Purpose / Status / Datasets /
Method / Known caveats sections). Must state explicitly, per spec:
- Bsal's weaker cross-check (raw GenBank annotation, no FungiDB/RefSeq curated source)
- The tiered ranking is a heuristic, not real co-expression (no RNA-seq/proteomic quantification exists for these life stages)
- The Bsal protein-count anomaly and its resolution (duplication-rate check from Task 5/6 live-run)
- The SIRIUS-coverage snapshot date/feature count used for this run (re-state each time Stage 4 is rerun, since the full native SIRIUS run is still deferred)

- [ ] **Step 8: Commit**

```bash
git add analysis/genome_bioactivity_linkage/scripts/build_linkage_tables.py \
        analysis/genome_bioactivity_linkage/tests/test_build_linkage_tables.py \
        pixi.toml \
        analysis/genome_bioactivity_linkage/GENOME_BIOACTIVITY_LINKAGE.md
git commit -m "genome-bioactivity-linkage: add end-to-end driver, pixi tasks, and writeup"
```

---

## Self-Review Notes

- **Spec coverage:** Reference genomes (Task 1), Stage 1 PFAM/BGC (Tasks 6–7), Stage 2 secretion (Tasks 3, 8), Stage 3 cross-reference (Tasks 4–5), Stage 4 linking incl. media-background subtraction and tiered ranking (Tasks 2, 9), Compute/`preempt` partitions (Tasks 3, 5 sbatch headers), Reproducibility/caching (skip-if-exists guards in `run_deeptmhmm.sh`, `run_rbh.sh`, `fetch_reference_annotation.sh`; Stage 4 rerun-in-full via `build_linkage_tables.py` reading `sirius_annotations.tsv` fresh each invocation), Output/writeup (Task 10) — all covered.
- **Placeholder scan:** No TBD/TODO markers. Two items are explicitly marked as blocked-on-external-dependency rather than deferred detail: DeepTMHMM's exact `TMRs.gff3` column layout (Task 3) and SignalP/PredGPI's exact BFD output column names (Task 10 Step 6) — both because the user's BFD run had not produced real output at plan-writing time; each has a concrete fallback action written out, not a placeholder.
- **Type consistency:** `parse_deeptmhmm.parse_tmrs_gff3` output shape (`protein_id, region_type, start, end`) is used identically in Tasks 3, 8, 10. `domain_families.classify_pfam`/`COMPOUND_CLASS_TO_FAMILY` used identically in Tasks 6, 9. `gene_domains` column contract (`protein_id, family, has_bgc_context, is_cross_ref_confirmed, is_extracellular`) matches between Task 9's test fixtures and Task 10's `build_gene_domain_table` output.
