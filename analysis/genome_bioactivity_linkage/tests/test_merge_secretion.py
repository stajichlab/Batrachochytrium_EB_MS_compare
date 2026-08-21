import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from merge_secretion import load_predgpi_gff3, load_signalp_gff3, predicted_extracellular  # noqa: E402

# Lines trimmed verbatim from a real BFD SignalP GFF3 output file:
# /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/function/
# signalp/9f/FF9F6419.signalp.gff3.gz
SIGNALP_GFF3_FIXTURE = """\
## gff-version 3
FF9F6419_000010-T1 FF9F6419_000010\tSignalP-6.0\tsignal_peptide\t1\t20\t0.9997437\t.\t.\t.
FF9F6419_000019-T1 FF9F6419_000019\tSignalP-6.0\tsignal_peptide\t1\t16\t0.53807586\t.\t.\t.
FF9F6419_000213-T1 FF9F6419_000213\tSignalP-6.0\tsignal_peptide\t1\t24\t0.60535\t.\t.\t.
"""

# Lines trimmed verbatim from a real BFD PredGPI GFF3 output file:
# /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/results/function/
# predgpi/9c/FE87067D.predgpi.gff3.gz
PREDGPI_GFF3_FIXTURE = """\
FE87067D_000001-T1\tPredGPI\tChain\t1\t431\t1.0\t.\t.\tevidence=ECO:0000256
FE87067D_000002-T1\tPredGPI\tChain\t1\t719\t1.0\t.\t.\tevidence=ECO:0000256
FE87067D_000017-T1\tPredGPI\tGPI-anchor\t408\t408\t1.0\t.\t.\tOntology_term=GO:0046658;evidence=ECO:0000256
"""


def test_load_signalp_gff3_parses_protein_id_and_cleavage_site(tmp_path):
    p = tmp_path / "signalp.gff3"
    p.write_text(SIGNALP_GFF3_FIXTURE)
    df = load_signalp_gff3(p)
    assert list(df.columns) == ["protein_id", "is_signal_peptide", "cleavage_site"]
    assert len(df) == 3
    # seqid column is "<protein_id> <protein_id>" -- must take only the first token.
    assert set(df["protein_id"]) == {
        "FF9F6419_000010-T1", "FF9F6419_000019-T1", "FF9F6419_000213-T1",
    }
    assert (df["is_signal_peptide"] == True).all()  # noqa: E712
    row = df[df["protein_id"] == "FF9F6419_000010-T1"].iloc[0]
    assert row["cleavage_site"] == 20


def test_load_signalp_gff3_gzip_transparent(tmp_path):
    import gzip

    p = tmp_path / "signalp.gff3.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(SIGNALP_GFF3_FIXTURE)
    df = load_signalp_gff3(p)
    assert len(df) == 3


def test_load_predgpi_gff3_derives_has_gpi_anchor_from_feature_type(tmp_path):
    p = tmp_path / "predgpi.gff3"
    p.write_text(PREDGPI_GFF3_FIXTURE)
    df = load_predgpi_gff3(p)
    assert list(df.columns) == ["protein_id", "has_gpi_anchor"]
    assert len(df) == 3
    by_id = df.set_index("protein_id")["has_gpi_anchor"].to_dict()
    assert bool(by_id["FE87067D_000001-T1"]) is False
    assert bool(by_id["FE87067D_000002-T1"]) is False
    assert bool(by_id["FE87067D_000017-T1"]) is True


def test_load_predgpi_gff3_gzip_transparent(tmp_path):
    import gzip

    p = tmp_path / "predgpi.gff3.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(PREDGPI_GFF3_FIXTURE)
    df = load_predgpi_gff3(p)
    assert len(df) == 3


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


def test_protein_missing_from_predgpi_is_not_falsely_extracellular():
    signalp = pd.DataFrame([{"protein_id": "protE", "is_signal_peptide": True, "cleavage_site": 25}])
    deeptmhmm = pd.DataFrame(
        [{"protein_id": "protE", "region_type": "signal", "start": 1, "end": 25}]
    )
    # protE is deliberately ABSENT from predgpi (not just has_gpi_anchor=False)
    predgpi = pd.DataFrame(columns=["protein_id", "has_gpi_anchor"])
    result = predicted_extracellular(signalp, deeptmhmm, predgpi)
    row = result[result["protein_id"] == "protE"].iloc[0]
    # Verify is_extracellular is proper bool dtype (numpy.bool_ or Python bool)
    assert isinstance(row["is_extracellular"], (bool, np.bool_))
    assert bool(row["is_extracellular"]) is True
