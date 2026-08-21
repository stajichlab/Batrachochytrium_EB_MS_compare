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
