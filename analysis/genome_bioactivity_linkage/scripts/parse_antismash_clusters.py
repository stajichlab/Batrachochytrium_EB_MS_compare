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
