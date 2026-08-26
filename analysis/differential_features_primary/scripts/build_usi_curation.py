#!/usr/bin/env python3
"""Curation scaffold for the top liquid-enriched bioactivity-flagged features.

Assembles, per species, the liq-vs-spore ENRICHED (log2FC_a_over_b > 0)
significant features with SIRIUS identity and (where present) the GNPS
spectral-library hit, ranks them by q_value then |log2FC|, and emits:

  - liq_enriched_curation/<species>_liq_enriched_top.tsv   (full ranked list)
  - liq_enriched_curation/<species>_liq_enriched_usi.html  (clickable USI grid)

Every row carries the live metabolomics-USI resolver link for its MS/MS and,
when a library hit exists, the node-vs-library mirror link, so each top
candidate can be visually verified against the spectrum pulled live from the
GNPS2 task (see data/raw/gnps2_e9838293_bagel/README_FOR_CLAUDE.md
"metabolomics-USI resolver" for the URL scheme).
"""
from __future__ import annotations

import html
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parents[3]
DIFF = REPO / "analysis" / "differential_features_primary"
OUT = DIFF / "liq_enriched_curation"
TASK = "e983829350de4bb39f278cbf22553247"
USI_BASE = f"mzspec:GNPS2:TASK-{TASK}-nf_output/feature_finding/aligned_features_filled.mgf:scan:"


def gene_uri(feature_id: int, suffix: str) -> str:
    # Trailing slash before the query string is required -- without it the
    # resolver 308-redirects to the http:// (not https://) version of the
    # same URL (confirmed 2026-08-26: `/json?usi1=...` -> 308 ->
    # `http://.../json/?usi1=...`), which browsers/img-tags/some HTTP
    # clients don't reliably follow (mixed-content blocking on an https
    # page, or a client that doesn't auto-follow), so links/images
    # silently fail to resolve to a specific spectrum.
    usi = f"{USI_BASE}{feature_id}"
    return f"https://metabolomics-usi.gnps2.org/{suffix}/?usi1={usi}"


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    diff = pd.read_csv(DIFF / "all_significant_features_summary.tsv", sep="\t")
    lib = pd.read_csv(
        REPO / "data" / "raw" / "gnps2_e9838293_bagel" / "nf_output" / "feature_library_search" / "merged_feature_library_search_results.tsv",
        sep="\t",
    )
    lib = lib.sort_values("cosine", ascending=False).drop_duplicates("query_scan")

    # all_significant_features_summary.tsv already carries the SIRIUS join
    # (sirius_structure_name/formula/npc_pathway/npc_class, annotation_origin,
    # bioactive flag) joined by differential-features-primary -- do NOT merge
    # sirius_annotations.tsv again here or the columns get _x/_y suffixed.
    rows = diff[diff["comparison"].str.contains(r"liq_\w+_vs_spore_\w+", regex=True) & (diff["log2FC_a_over_b"] > 0)].copy()
    rows["abs_log2fc"] = rows["log2FC_a_over_b"].abs()
    libs = lib.rename(columns={"query_scan": "row_id"})
    rows = rows.merge(libs[["row_id", "NAME", "FORMULA", "cosine", "matched_peaks", "SPECTRUMID"]], on="row_id", how="left")
    rows = rows.sort_values(["q_value", "abs_log2fc"], ascending=[True, False])
    rows = rows.drop_duplicates("row_id")

    for species in ["dendrobatidis", "salamandrivorans"]:
        df = rows[rows["species"] == species].drop_duplicates(subset=["row_id"]).copy()
        df["abs_log2fc"] = df["log2FC_a_over_b"].abs()
        df = df.sort_values(["q_value", "abs_log2fc"], ascending=[True, False])
        cols = [
            "row_id", "mz", "rt", "comparison", "log2FC_a_over_b", "q_value",
            "bioactive", "annotation_origin",
            "sirius_structure_name", "sirius_formula", "sirius_npc_pathway", "sirius_npc_class",
            "NAME", "cosine", "matched_peaks",
        ]
        cols = [c for c in cols if c in df.columns]
        df = df[cols]
        df.to_csv(OUT / f"{species}_liq_enriched_top.tsv", sep="\t", index=False)
        write_html(df, species)

        # aggregate stats for the writeup
        n = len(df)
        n_struc = int(df["sirius_structure_name"].notna().sum())
        n_lib = int(df["NAME"].notna().sum())
        n_nprs = int(df["sirius_npc_pathway"].isin(["Amino acids and Peptides", "Alkaloids", "Polyketides", "Terpenoids"]).sum())
        print(f"{species}: {n} unique liq-enriched features, {n_struc} SIRIUS-structure, {n_lib} GNPS-library-hit, {n_nprs} in mapped compound classes")


def write_html(df: pd.DataFrame, species: str) -> None:
    n = min(len(df), 100)
    top = df.head(n)
    rows_html = []
    for _, r in top.iterrows():
        fid = int(r["row_id"])
        img = gene_uri(fid, "png")
        spec = gene_uri(fid, "spectrum")
        js = gene_uri(fid, "json")
        name = r.get("sirius_structure_name") or r.get("NAME") or ""
        name = html.escape(str(name)) if name else "&mdash;"
        lib_usi = ""
        if pd.notna(r.get("SPECTRUMID")):
            lid = r["SPECTRUMID"]
            mirror = f"https://metabolomics-usi.gnps2.org/png/mirror/?usi1={USI_BASE}{fid}&usi2=mzspec:GNPS:GNPS-LIBRARY:accession:{lid}"
            lib_usi = f'<a href="{mirror}" target="_blank">&#8617;vs-library ({lid})</a>'
        rows_html.append(
            "<tr>"
            f'<td><a href="{spec}" target="_blank">{fid}</a></td>'
            f"<td>{r['mz']:.4f}</td>"
            f"<td>{r['rt']:.2f}</td>"
            f"<td>{r['log2FC_a_over_b']:.1f}</td>"
            f"<td>{r['q_value']:.2e}</td>"
            f"<td>{name}</td>"
            f"<td>{r.get('sirius_formula','') or ''}</td>"
            f'<td><a href="{img}" target="_blank"><img src="{img}" style="max-width:200px;max-height:90px"></a></td>'
            f'<td><a href="{spec}" target="_blank">view</a> / <a href="{js}" target="_blank">json</a> {lib_usi}</td>'
            "</tr>"
        )
    html_doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>{species} liq-enriched curation</title></head>
<body><h1>{species} &mdash; top liq-vs-spore-enriched features with live USI spectra</h1>
<p>Top {n} of {len(df)} unique liq-enriched significant features, ranked by q_value then |log2FC|.</p>
<table border="1" cellspacing="0" cellpadding="4">
<tr><th>feature</th><th>m/z</th><th>RT(min)</th><th>log2FC</th><th>q</th><th>structure</th><th>formula</th><th>spectrum</th><th>actions</th></tr>
{''.join(rows_html)}
</table></body></html>"""
    (OUT / f"{species}_liq_enriched_usi.html").write_text(html_doc)


if __name__ == "__main__":
    main()
