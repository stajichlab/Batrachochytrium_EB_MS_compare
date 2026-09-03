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

MEDIA-BLANK ORDERING (2026-09-02): "liq-enriched" alone is NOT evidence of
secretion. Both media are peptide-rich broths (Bd 1% tryptone = casein
digest, Bsal 50% TGHL), so medium-derived peptides are abundant in `liq` and
absent from the washed `spore` pellet and therefore score as maximally
liq-enriched. Only ~9% (Bd) / ~20% (Bsal) of the rows here exceed their own
C_liq blank. Rows are consequently ORDERED by `passes_media_blank` first, so
the HTML grid -- the actual MS2 shortlist -- shows blank-clearing features.
No rows are dropped (see MYCELIUM.md: do not subset without confirmation);
`passes_media_blank` is a column on the TSV so the full list stays available.
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


def _precursor_concordance(bagel: pd.DataFrame, tol: float = 0.01) -> dict[int, bool]:
    """row_id -> does the MGF precursor agree with the feature table's m/z?

    Added 2026-09-02 after an audit found annotations built on the wrong
    precursor. Parsing all MGF blocks and comparing PEPMASS to `row m/z` for
    the 6,453 has_ms2 features: 6,199 agree within 0.01, but 76 differ by
    0.01-0.3, 57 by 0.3-0.7 (half-integer -- the signature of a 2+ ion written
    into the table as charge=1/M+0/default-adduct, which the artifact filter
    therefore cannot see) and 105 by ~1.0 (isotope spacing). The discordance is
    essentially confined to the 489 SOURCE_FEATURE_ID=-1 blocks.

    SIRIUS was handed those wrong precursors verbatim, which is where the
    chemically impossible shortlist formulas come from (C10H5Cl9,
    C16H21Br4N3O2, C21H19Br2IN6O6 -- in a fungal culture in tryptone). Those
    are not annotations, they are the formula finder absorbing a mass error.
    """
    mgf = (REPO / "data" / "raw" / "gnps2_e9838293_bagel" / "nf_output"
           / "feature_finding" / "aligned_features_filled.mgf")
    pepmass: dict[int, float] = {}
    cur: dict[str, str] = {}
    with open(mgf) as fh:
        for line in fh:
            t = line.strip()
            if t == "BEGIN IONS":
                cur = {}
            elif t == "END IONS":
                if "SCANS" in cur and "PEPMASS" in cur:
                    try:
                        pepmass[int(cur["SCANS"])] = float(cur["PEPMASS"].split()[0])
                    except ValueError:
                        pass
                cur = {}
            elif "=" in t and t.split("=")[0].isupper():
                k, v = t.split("=", 1)
                cur[k] = v
    out: dict[int, bool] = {}
    for rid, mz in zip(bagel["row_id"], bagel["row_mz"]):
        pm = pepmass.get(int(rid))
        # No MGF block -> cannot verify; has_ms2 gates those out anyway.
        out[int(rid)] = pm is not None and abs(pm - mz) <= tol
    return out


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
    # has_ms2 is NOT in linked_data; pull it from the bagel table. Without it
    # the "MS2 verification shortlist" is largely unverifiable -- only 6,453 of
    # 38,547 features have an acquired MS2 spectrum, and the USI resolver
    # renders the GAP-FILLED mgf regardless, so a no-MS2 row still returns a
    # picture. 73% of the pre-2026-09-02 top-100 grid had no MS2 at all.
    bagel = pd.read_csv(
        REPO / "data" / "raw" / "gnps2_e9838293_bagel" / "nf_output" / "feature_finding"
        / "feature_finding_results" / "aligned_features.csv",
        low_memory=False, usecols=["row ID", "row m/z", "has_ms2"],
    ).rename(columns={"row ID": "row_id", "row m/z": "row_mz"})
    bagel["mz_concordant"] = bagel["row_id"].map(_precursor_concordance(bagel))

    rows = diff[diff["comparison"].str.contains(r"liq_\w+_vs_spore_\w+", regex=True) & (diff["log2FC_a_over_b"] > 0)].copy()
    rows["abs_log2fc"] = rows["log2FC_a_over_b"].abs()
    libs = lib.rename(columns={"query_scan": "row_id"})
    rows = rows.merge(libs[["row_id", "NAME", "FORMULA", "cosine", "matched_peaks", "SPECTRUMID"]], on="row_id", how="left")
    rows = rows.merge(bagel[["row_id", "has_ms2", "mz_concordant"]], on="row_id", how="left")
    rows["has_ms2"] = rows["has_ms2"].fillna(False).astype(bool)
    rows["mz_concordant"] = rows["mz_concordant"].fillna(False).astype(bool)
    rows = rows.sort_values(["q_value", "abs_log2fc"], ascending=[True, False])
    rows = rows.drop_duplicates("row_id")

    for species in ["dendrobatidis", "salamandrivorans"]:
        df = rows[rows["species"] == species].drop_duplicates(subset=["row_id"]).copy()
        df["abs_log2fc"] = df["log2FC_a_over_b"].abs()
        if "passes_media_blank" not in df.columns:
            raise SystemExit(
                "all_significant_features_summary.tsv has no passes_media_blank column -- "
                "regenerate it with `pixi run differential-features-primary && "
                "pixi run feature-tables-primary` first (column added 2026-09-02)."
            )
        df["passes_media_blank"] = df["passes_media_blank"].fillna(False).astype(bool)
        # Blank-clearing features first; the HTML grid is the MS2 shortlist.
        # Shortlist order: must clear the media blank AND have an acquired MS2
        # before q-value matters. Ranking within the top by |log2FC| would be
        # meaningless anyway -- q is tied at the attainable floor across the
        # whole head of the list, so |log2FC| silently decided the order, and
        # log2FC is pseudocount-bounded for on/off features.
        # mz_concordant added 2026-09-02: an annotation built on the wrong
        # precursor is not a weaker annotation, it is a different molecule.
        df["shortlist_ready"] = (
            df["passes_media_blank"] & df["has_ms2"] & df["mz_concordant"]
        )
        df = df.sort_values(
            ["shortlist_ready", "passes_media_blank", "has_ms2", "mz_concordant",
             "q_value", "prevalence_diff"],
            ascending=[False, False, False, False, True, False],
        )
        cols = [
            "row_id", "mz", "rt", "comparison", "log2FC_a_over_b", "q_value",
            "shortlist_ready", "passes_media_blank", "has_ms2", "mz_concordant",
            "prevalence_diff",
            "is_liq_enriched", "bioactive", "annotation_origin",
            "sirius_structure_confidence",
            "sirius_structure_name", "sirius_formula", "sirius_npc_pathway", "sirius_npc_class",
            "NAME", "cosine", "matched_peaks",
        ]
        cols = [c for c in cols if c in df.columns]
        df = df[cols]
        df.to_csv(OUT / f"{species}_liq_enriched_top.tsv", sep="\t", index=False)
        write_html(df, species)

        # aggregate stats for the writeup
        n = len(df)
        n_blank = int(df["passes_media_blank"].sum())
        n_ready = int(df["shortlist_ready"].sum())
        n_struc = int(df["sirius_structure_name"].notna().sum())
        n_lib = int(df["NAME"].notna().sum())
        n_nprs = int(df["sirius_npc_pathway"].isin(["Amino acids and Peptides", "Alkaloids", "Polyketides", "Terpenoids"]).sum())
        print(
            f"{species}: {n} unique liq-enriched features "
            f"({n_blank}, {n_blank / n:.1%}, clear the C_liq media blank; "
            f"{n_ready} ALSO have acquired MS2 = shortlist-ready), "
            f"{n_struc} SIRIUS-structure, {n_lib} GNPS-library-hit, {n_nprs} in mapped compound classes"
        )


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
        blank_ok = bool(r.get("passes_media_blank", False))
        # Amber-tint any row that does not clear its own media blank, so a
        # medium-derived peptide can never be mistaken for a secreted product.
        tr_style = "" if blank_ok else ' style="background:#FFF4E5"'
        blank_cell = "&#10003;" if blank_ok else "&#10007; media"
        rows_html.append(
            f"<tr{tr_style}>"
            f'<td><a href="{spec}" target="_blank">{fid}</a></td>'
            f"<td>{r['mz']:.4f}</td>"
            f"<td>{r['rt']:.2f}</td>"
            f"<td>{r['log2FC_a_over_b']:.1f}</td>"
            f"<td>{r['q_value']:.2e}</td>"
            f"<td>{blank_cell}</td>"
            f"<td>{name}</td>"
            f"<td>{r.get('sirius_formula','') or ''}</td>"
            f'<td><a href="{img}" target="_blank"><img src="{img}" style="max-width:200px;max-height:90px"></a></td>'
            f'<td><a href="{spec}" target="_blank">view</a> / <a href="{js}" target="_blank">json</a> {lib_usi}</td>'
            "</tr>"
        )
    n_blank_total = int(df["passes_media_blank"].sum())
    n_blank_shown = int(top["passes_media_blank"].sum())
    html_doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>{species} liq-enriched curation</title></head>
<body><h1>{species} &mdash; top liq-vs-spore-enriched features with live USI spectra</h1>
<p>Top {n} of {len(df)} unique liq-enriched significant features, ranked by
media-blank status, then q_value, then |log2FC|.
<b>{n_blank_shown} of the {n} shown</b> (and {n_blank_total} of {len(df)} overall)
clear the C_liq media blank at &ge;2&times;.</p>
<p style="background:#FFF4E5;padding:6px;border:1px solid #E8C48A">
<b>Read the blank column first.</b> Both media are peptide-rich broths
(Bd 1% tryptone = casein digest, Bsal 50% TGHL). Medium-derived peptides are
abundant in <code>liq</code> and absent from the washed <code>spore</code>
pellet, so they score as maximally &quot;liq-enriched&quot; without being
secreted at all. Amber rows (&#10007; media) do <em>not</em> exceed their own
media blank and are not evidence of secretion.</p>
<table border="1" cellspacing="0" cellpadding="4">
<tr><th>feature</th><th>m/z</th><th>RT(min)</th><th>log2FC</th><th>q</th><th>&gt;2&times; blank</th><th>structure</th><th>formula</th><th>spectrum</th><th>actions</th></tr>
{''.join(rows_html)}
</table></body></html>"""
    (OUT / f"{species}_liq_enriched_usi.html").write_text(html_doc)


if __name__ == "__main__":
    main()
