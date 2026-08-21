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
