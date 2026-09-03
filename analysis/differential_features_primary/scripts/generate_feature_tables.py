#!/usr/bin/env python3
"""
Feature-table HTML for investigating the primary-tier differential results
(analysis/differential_features_primary/). Port of the Rhodotorula
generate_compound_table_html.py strategy to this project's schema.

For each of the 8 primary contrasts it turns the significant (q < 0.05)
feature rows into a single self-contained, offline-viewable
`compound_summary.html` (data embedded as JSON, sort/filter in plain JS, no
CDN). The rollup is written as one TSV plus one per-species HTML
(`all_significant_features_summary_<species>.html`) — the combined HTML would
be ~125 MB, over GitHub's 100 MB per-file hard limit, so it is split on the
`species` column so each chunk stays pushable.

Design decisions (inherited from the Rhodotorula generator):
  - Numeric columns sort numerically (q-value in sci notation / negative
    log2FC would sort wrong lexicographically); the raw value is embedded in
    the row object and only the rendered <td> is formatted.
  - Identity is a shape+color glyph in the leftmost data column (Okabe-Ito):
    filled diamond = SIRIUS structure, hollow circle = SIRIUS formula only,
    dash = unidentified.
  - ~10 columns in the main table; the other ~25 (SIRIUS/CANOPUS details,
    transfer match stats) are in a per-row expand panel on click.
  - Filters, in priority order: click-header sort, numeric thresholds
    (q-value <=, |log2FC| >=), identity-source chips, then checkboxes for
    the curated flags (bioactive, secreted candidate), then free-text search.
  - Plain <table> + DocumentFragment rendering with debounced search input.

Identity source per row is derived from the SIRIUS columns:
  sirius_structure_name non-null -> sirius_structure
  else sirius_formula non-null   -> sirius_formula_only
  else                           -> unidentified

Usage:
    python3 scripts/generate_feature_tables.py [--significant significant_annotated.tsv]

Writes:
  <comparison>/compound_summary.tsv        per-comparison significant table
  <comparison>/compound_summary.html       sortable/filterable view
  all_significant_features_summary.tsv     rollup across all 8 contrasts
  all_significant_features_summary_<species>.html
                                           rollup sortable/filterable view, chunked
                                           per species (dendrobatidis / salamandrivorans)
                                           so each file stays under GitHub's 100MB
                                           per-file hard limit.
  feature_tables_index.html                navigation hub linking the above
"""
import argparse
import html
import json
import math
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
OUT_ROOT = ROOT
SIG_TSV = OUT_ROOT / "significant_annotated.tsv"

# (json_key, header_label, kind) -- kind in {"int","float","sci","bool","text","glyph"}
PRIMARY_COLS = [
    ("best_identity_source", "ID", "glyph"),
    ("row_id", "row ID", "int"),
    ("mz", "m/z", "float"),
    ("rt", "RT (min)", "float"),
    ("log2FC_a_over_b", "log2FC (A/B)", "float"),
    ("q_value", "q-value", "sci"),
    ("direction", "direction", "text"),
    ("bioactive", "bioact", "bool"),
    ("is_secreted_candidate", "secreted", "bool"),
    ("best_identity", "best identity", "text"),
]
# `comparison` is only present in the rollup table; added there.
ROLLUP_COLS = [("comparison", "comparison", "text")] + PRIMARY_COLS

SECONDARY_LABELS = {
    "row_id": "feature row ID",
    "mz": "feature m/z",
    "rt": "RT (min)",
    "median_a": "median area (A)",
    "median_b": "median area (B)",
    "log2FC_a_over_b": "log2(median A / median B)",
    "U_stat": "Mann-Whitney U",
    "p_value": "Mann-Whitney p-value",
    "q_value": "BH-FDR q-value",
    "comparison": "contrast",
    "species": "species",
    "family": "contrast family",
    "group_a": "group A (numerator)",
    "group_b": "group B (denominator)",
    "direction": "direction (up in which group)",
    "liq_over_spore_log2fc": "liq-vs-spore log2FC hint (stage-confounded)",
    "is_secreted_candidate": "secreted/enriched-in-liquid candidate",
    "row ID": "SIRIUS row ID",
    "sirius_formula": "SIRIUS formula",
    "sirius_adduct": "SIRIUS adduct",
    "sirius_structure_name": "SIRIUS structure name",
    "sirius_structure_smiles": "SIRIUS structure SMILES",
    "sirius_structure_confidence": "SIRIUS structure confidence",
    "sirius_npc_pathway": "NPC pathway",
    "sirius_npc_class": "NPC class",
    "sirius_classyfire_class": "ClassyFire class",
    "sirius_source_feature_id": "SIRIUS source feature id",
    "sirius_source_run": "SIRIUS source run",
    "annotation_origin": "annotation origin (transferred/native)",
    "n_candidates": "SIRIUS candidate count",
    "source_mz": "source m/z",
    "source_rt": "source RT (min)",
    "feature_mz": "feature m/z",
    "feature_rt": "feature RT (min)",
    "ppm_error": "m/z ppm error",
    "rt_delta_min": "RT delta (min)",
    "ms2_cosine": "MS2 cosine",
    "match_class": "match class",
    "match_status": "match status",
    "n_sirius_hits": "SIRIUS hit count",
    "sirius_hit_ids": "SIRIUS hit IDs",
    "sirius_hit_formulas": "SIRIUS hit formulas",
    "merged_conflict": "merged annotation conflict",
    "bioactive": "bioactivity keyword flag",
}

GLYPH_BY_SOURCE = {
    "sirius_structure": {"glyph": "◇", "color": "#E69F00", "label": "SIRIUS structure"},
    "sirius_formula_only": {"glyph": "○", "color": "#999999", "label": "SIRIUS formula only"},
    "unidentified": {"glyph": "—", "color": "#BBBBBB", "label": "unidentified"},
}


def clean_value(v):
    if v is None:
        return None
    if isinstance(v, float) and (math.isnan(v) or math.isinf(v)):
        return None
    if isinstance(v, bool):
        return v
    return v


# --- payload compression -------------------------------------------------
# The embedded payload was an array-of-objects, which repeats every column
# NAME on every row: at 46 columns x 54,770 rows that is ~31 MB of key text
# in the salamandrivorans rollup alone, over half the file. Two lossless-for-
# display transforms shrink it with no runtime dependency and no change to
# the consuming JS (the decoder rebuilds the identical array-of-objects):
#
#   1. columnar layout    -- one array per column; each key appears once.
#   2. dictionary coding  -- low-cardinality columns (comparison, direction,
#                            family, species, group_a/b, npc_pathway, ...)
#                            become {values, integer indices}. `comparison`
#                            alone has 4 distinct values across 54,770 rows.
#
# Floats are rounded to FLOAT_SIG_DIGITS significant digits. This is a VIEW
# artifact -- the table renders at .4f (m/z), .2f (RT) and .2e (q) -- and the
# TSVs beside it remain the full-precision record.
FLOAT_SIG_DIGITS = 7
# Dictionary-encode only when the dictionary is a real win: the column must
# repeat substantially, or the indices cost more than the values they replace.
DICT_MAX_RATIO = 0.5


def _round_sig(v: float, sig: int = FLOAT_SIG_DIGITS):
    if v == 0 or not math.isfinite(v):
        return v
    r = round(v, sig - 1 - math.floor(math.log10(abs(v))))
    # Collapse 3.0 -> 3 so ints don't carry a redundant ".0" in the JSON.
    return int(r) if r == int(r) and abs(r) < 1e15 else r


def encode_columnar(df: pd.DataFrame, columns: list[str]) -> dict:
    """Column-oriented, dictionary-coded payload. Decoded by DECODER_JS."""
    out: dict[str, dict] = {}
    n = len(df)
    for key in columns:
        vals = [clean_value(v) for v in df[key]]
        vals = [_round_sig(v) if isinstance(v, float) else v for v in vals]
        # `None` participates in the dictionary like any other value.
        uniq = list(dict.fromkeys(map(_hashable, vals)))
        if n and len(uniq) <= max(1, int(n * DICT_MAX_RATIO)):
            index = {u: i for i, u in enumerate(uniq)}
            out[key] = {"t": "d", "v": uniq, "i": [index[_hashable(v)] for v in vals]}
        else:
            out[key] = {"t": "r", "v": vals}
    return {"n": n, "cols": out}


def _hashable(v):
    # Every value we emit is already a JSON scalar; this guards against a
    # stray list/dict silently becoming unhashable during dict-coding.
    return v if isinstance(v, (str, int, float, bool, type(None))) else str(v)


DECODER_JS = """
const DATA = (() => {
  const n = PAYLOAD.n, cols = PAYLOAD.cols, out = new Array(n);
  for (let r = 0; r < n; r++) out[r] = {};
  for (const key in cols) {
    const c = cols[key];
    if (c.t === 'd') {
      const v = c.v, ix = c.i;
      for (let r = 0; r < n; r++) out[r][key] = v[ix[r]];
    } else {
      const v = c.v;
      for (let r = 0; r < n; r++) out[r][key] = v[r];
    }
  }
  return out;
})();
"""


def identity_source(row):
    if pd.notna(row.get("sirius_structure_name")):
        return "sirius_structure"
    if pd.notna(row.get("sirius_formula")):
        return "sirius_formula_only"
    return "unidentified"


def best_identity(row):
    for col in ("sirius_structure_name", "sirius_formula"):
        v = row.get(col)
        if pd.notna(v) and str(v).strip():
            return str(v)
    return None


def build_html(df: pd.DataFrame, title: str, is_rollup: bool = False) -> str:
    df = df.copy()
    df["best_identity_source"] = df.apply(identity_source, axis=1)
    df["best_identity"] = df.apply(best_identity, axis=1)

    cols = ROLLUP_COLS if is_rollup else PRIMARY_COLS
    cols_present = [c for c in cols if c[0] in df.columns]
    all_json_cols = list(df.columns)
    secondary_cols = [c for c in all_json_cols if c not in {c0 for c0, *_ in PRIMARY_COLS}]

    payload = encode_columnar(df, all_json_cols)

    n_total = len(df)

    def _count_true(col: str) -> int:
        return int(df[col].fillna(False).astype(bool).sum()) if col in df.columns else 0

    n_identified = (
        int((~df["best_identity_source"].isin([None, "unidentified"]) & df["best_identity_source"].notna()).sum())
        if "best_identity_source" in df.columns else 0
    )
    n_bioactive = _count_true("bioactive")
    n_secreted = _count_true("is_secreted_candidate")

    primary_cols_json = [{"key": k, "label": lbl, "kind": kind} for k, lbl, kind in cols_present]
    secondary_labels_json = {k: SECONDARY_LABELS.get(k, k) for k in secondary_cols}

    return TEMPLATE.format(
        title=html.escape(title),
        n_total=n_total,
        n_identified=n_identified,
        n_bioactive=n_bioactive,
        n_secreted=n_secreted,
        data_json=json.dumps(payload, separators=(",", ":")),
        decoder_js=DECODER_JS,
        primary_cols_json=json.dumps(primary_cols_json),
        secondary_cols_json=json.dumps(secondary_cols),
        secondary_labels_json=json.dumps(secondary_labels_json),
        glyphs_json=json.dumps(GLYPH_BY_SOURCE),
    )


TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>
  :root {{
    --bg: #ffffff; --surface: #f7f7f8; --border: #e2e2e6;
    --ink: #1a1a1e; --ink-muted: #6b6b74; --accent: #0072B2;
    --row-hover: #f0f4f8; --flag-bg: #eef8f2;
  }}
  * {{ box-sizing: border-box; }}
  body {{
    margin: 0; padding: 24px; background: var(--bg); color: var(--ink);
    font: 14px/1.4 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
  }}
  h1 {{ font-size: 18px; margin: 0 0 4px; }}
  .subtitle {{ color: var(--ink-muted); margin: 0 0 16px; font-size: 13px; }}
  .toolbar {{
    display: flex; flex-wrap: wrap; gap: 10px; align-items: center;
    padding: 12px; background: var(--surface); border: 1px solid var(--border);
    border-radius: 8px; margin-bottom: 12px; position: sticky; top: 0; z-index: 2;
  }}
  .toolbar label {{ font-size: 12px; color: var(--ink-muted); display: flex; align-items: center; gap: 6px; }}
  .toolbar input[type="text"], .toolbar input[type="number"] {{
    font: inherit; padding: 5px 8px; border: 1px solid var(--border); border-radius: 5px; width: 100px;
  }}
  .toolbar input[type="text"] {{ width: 220px; }}
  .chip {{
    font: inherit; font-size: 12px; padding: 5px 10px; border-radius: 999px;
    border: 1px solid var(--border); background: #fff; cursor: pointer; color: var(--ink);
  }}
  .chip[aria-pressed="true"] {{ background: var(--accent); color: #fff; border-color: var(--accent); }}
  .chip .swatch {{ display: inline-block; width: 10px; text-align: center; margin-right: 4px; }}
  #count {{ margin-left: auto; font-size: 12px; color: var(--ink-muted); }}
  table {{ border-collapse: collapse; width: 100%; font-variant-numeric: tabular-nums; }}
  thead th {{
    position: sticky; top: 58px; background: var(--surface); text-align: left;
    padding: 8px 10px; border-bottom: 2px solid var(--border); cursor: pointer;
    white-space: nowrap; user-select: none; font-size: 12px; color: var(--ink-muted);
  }}
  thead th:hover {{ color: var(--ink); }}
  thead th.sorted {{ color: var(--accent); }}
  thead th .arrow {{ font-size: 10px; margin-left: 3px; }}
  tbody td {{ padding: 7px 10px; border-bottom: 1px solid var(--border); white-space: nowrap; }}
  tbody tr.data-row {{ cursor: pointer; }}
  tbody tr.data-row:hover {{ background: var(--row-hover); }}
  tbody tr.data-row.flagged {{ background: var(--flag-bg); }}
  tbody tr.data-row.flagged:hover {{ background: #e2f2ea; }}
  td.identity {{ text-align: center; font-size: 15px; }}
  td.name-cell {{ max-width: 380px; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
  td.flag-cell {{ text-align: center; }}
  tr.detail-row {{ display: none; }}
  tr.detail-row.open {{ display: table-row; }}
  tr.detail-row td {{ background: var(--surface); white-space: normal; padding: 12px 16px; }}
  .detail-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(260px, 1fr)); gap: 6px 20px; }}
  .detail-item dt {{ font-size: 11px; color: var(--ink-muted); margin: 0; }}
  .detail-item dd {{ margin: 0 0 6px; font-size: 13px; word-break: break-word; }}
  .legend {{ display: flex; gap: 14px; flex-wrap: wrap; font-size: 12px; color: var(--ink-muted); margin: 0 0 12px; }}
  .legend span.g {{ margin-right: 3px; }}
  @media (prefers-color-scheme: dark) {{
    :root {{
      --bg: #16161a; --surface: #1e1e23; --border: #303038;
      --ink: #eaeaef; --ink-muted: #9a9aa5; --row-hover: #26262c; --flag-bg: #17251d;
    }}
    .chip {{ background: #222228; }}
    tr.detail-row.flagged:hover {{ background: #1c2f24; }}
  }}
</style>
</head>
<body>
<h1>{title}</h1>
<p class="subtitle">{n_total} significant features &middot; {n_identified} SIRIUS-annotated &middot; {n_bioactive} bioactivity-flagged &middot; {n_secreted} secreted-candidates</p>
<p class="legend" id="legend"></p>
<div class="toolbar">
  <label>Search <input type="text" id="search" placeholder="identity, formula, adduct..."></label>
  <label>q-value &le; <input type="number" id="qmax" step="0.001" placeholder="0.05"></label>
  <label>|log2FC| &ge; <input type="number" id="fcmin" step="0.1" placeholder="0"></label>
  <label><input type="checkbox" id="bioactOnly"> bioact only</label>
  <label><input type="checkbox" id="secretedOnly"> secreted only</label>
  <span id="sourceChips"></span>
  <span id="count"></span>
</div>
<table>
  <thead><tr id="headRow"></tr></thead>
  <tbody id="body"></tbody>
</table>

<script>
const PAYLOAD = {data_json};
{decoder_js}
const PRIMARY_COLS = {primary_cols_json};
const SECONDARY_COLS = {secondary_cols_json};
const SECONDARY_LABELS = {secondary_labels_json};
const GLYPHS = {glyphs_json};

let sortKey = "q_value", sortDir = 1;
let openRow = null;

function fmt(kind, v) {{
  if (v === null || v === undefined) return "–";
  if (kind === "sci") return Number(v).toExponential(2);
  if (kind === "bool") return v ? "✓" : "";
  return String(v);
}}
function fmtFloat(v, digits) {{
  if (v === null || v === undefined) return "–";
  return Number(v).toFixed(digits);
}}

function glyphCell(source) {{
  const g = GLYPHS[source] || GLYPHS["unidentified"];
  return `<span title="${{g.label}}" style="color:${{g.color}}">${{g.glyph}}</span>`;
}}

function buildLegend() {{
  const el = document.getElementById("legend");
  el.innerHTML = Object.values(GLYPHS).map(g =>
    `<span><span class="g" style="color:${{g.color}}">${{g.glyph}}</span>${{g.label}}</span>`
  ).join("");
}}

function buildHead() {{
  const tr = document.getElementById("headRow");
  tr.innerHTML = "";
  PRIMARY_COLS.forEach(col => {{
    const th = document.createElement("th");
    th.textContent = col.label;
    th.dataset.key = col.key;
    th.dataset.kind = col.kind;
    const arrow = document.createElement("span");
    arrow.className = "arrow";
    th.appendChild(arrow);
    th.addEventListener("click", () => {{
      if (sortKey === col.key) sortDir *= -1; else {{ sortKey = col.key; sortDir = col.kind === "text" || col.kind === "glyph" ? 1 : -1; }}
      render();
    }});
    tr.appendChild(th);
  }});
  const th = document.createElement("th");
  tr.appendChild(th); // spacer for detail toggle
}}

let sourceFilter = new Set(Object.keys(GLYPHS));
function buildSourceChips() {{
  const el = document.getElementById("sourceChips");
  el.innerHTML = "";
  Object.entries(GLYPHS).forEach(([key, g]) => {{
    const btn = document.createElement("button");
    btn.className = "chip";
    btn.setAttribute("aria-pressed", "true");
    btn.innerHTML = `<span class="swatch" style="color:${{g.color}}">${{g.glyph}}</span>${{g.label}}`;
    btn.addEventListener("click", () => {{
      if (sourceFilter.has(key)) {{ sourceFilter.delete(key); btn.setAttribute("aria-pressed", "false"); }}
      else {{ sourceFilter.add(key); btn.setAttribute("aria-pressed", "true"); }}
      render();
    }});
    el.appendChild(btn);
  }});
}}

function currentFilters() {{
  return {{
    search: document.getElementById("search").value.trim().toLowerCase(),
    qmax: parseFloat(document.getElementById("qmax").value),
    fcmin: parseFloat(document.getElementById("fcmin").value),
    bioactOnly: document.getElementById("bioactOnly").checked,
    secretedOnly: document.getElementById("secretedOnly").checked,
  }};
}}

function matches(row, f) {{
  if (!sourceFilter.has(row.best_identity_source || "unidentified")) return false;
  if (f.bioactOnly && row.bioactive !== true) return false;
  if (f.secretedOnly && row.is_secreted_candidate !== true) return false;
  if (!isNaN(f.qmax) && !(row.q_value <= f.qmax)) return false;
  if (!isNaN(f.fcmin) && !(Math.abs(row.log2FC_a_over_b) >= f.fcmin)) return false;
  if (f.search) {{
    const hay = [row.row_id, row.best_identity, row.sirius_structure_name,
                 row.sirius_formula, row.sirius_adduct, row.comparison].filter(Boolean).join(" ").toLowerCase();
    if (!hay.includes(f.search)) return false;
  }}
  return true;
}}

function detailPanel(row) {{
  const items = SECONDARY_COLS.map(k => {{
    const v = row[k];
    if (v === null || v === undefined || v === "") return "";
    return `<div class="detail-item"><dt>${{SECONDARY_LABELS[k] || k}}</dt><dd>${{String(v)}}</dd></div>`;
  }}).join("");
  return `<dl class="detail-grid">${{items}}</dl>`;
}}

let debounceTimer;
function scheduleRender() {{
  clearTimeout(debounceTimer);
  debounceTimer = setTimeout(render, 150);
}}

function render() {{
  const f = currentFilters();
  let rows = DATA.filter(r => matches(r, f));

  document.querySelectorAll("#headRow th[data-key]").forEach(th => {{
    th.classList.toggle("sorted", th.dataset.key === sortKey);
    const arrow = th.querySelector(".arrow");
    arrow.textContent = th.dataset.key === sortKey ? (sortDir === 1 ? "▲" : "▼") : "";
  }});

  const kind = (PRIMARY_COLS.find(c => c.key === sortKey) || {{}}).kind;
  rows.sort((a, b) => {{
    let av = a[sortKey], bv = b[sortKey];
    if (kind === "text" || kind === "glyph") {{
      av = (av || "").toString(); bv = (bv || "").toString();
      return sortDir * av.localeCompare(bv);
    }}
    av = av === null || av === undefined ? -Infinity : Number(av);
    bv = bv === null || bv === undefined ? -Infinity : Number(bv);
    return sortDir * (av - bv);
  }});

  document.getElementById("count").textContent = `${{rows.length}} / ${{DATA.length}} rows`;

  const tbody = document.getElementById("body");
  tbody.innerHTML = "";
  const frag = document.createDocumentFragment();

  rows.forEach((row, i) => {{
    const tr = document.createElement("tr");
    tr.className = "data-row" + (row.bioactive ? " flagged" : "");
    PRIMARY_COLS.forEach(col => {{
      const td = document.createElement("td");
      if (col.key === "best_identity_source") {{
        td.className = "identity";
        td.innerHTML = glyphCell(row.best_identity_source);
      }} else if (col.key === "best_identity") {{
        td.className = "name-cell";
        td.textContent = row.best_identity || "–";
      }} else if (col.kind === "bool") {{
        td.className = "flag-cell";
        td.textContent = row[col.key] ? "✓" : "";
      }} else if (col.kind === "float") {{
        td.textContent = fmtFloat(row[col.key], col.key === "rt" ? 2 : 4);
      }} else {{
        td.textContent = fmt(col.kind, row[col.key]);
      }}
      tr.appendChild(td);
    }});
    const toggleTd = document.createElement("td");
    toggleTd.textContent = "›";
    toggleTd.style.color = "var(--ink-muted)";
    tr.appendChild(toggleTd);

    const detailTr = document.createElement("tr");
    detailTr.className = "detail-row";
    const detailTd = document.createElement("td");
    detailTd.colSpan = PRIMARY_COLS.length + 1;
    detailTd.innerHTML = detailPanel(row);
    detailTr.appendChild(detailTd);

    tr.addEventListener("click", () => {{
      const wasOpen = detailTr.classList.contains("open");
      document.querySelectorAll("tr.detail-row.open").forEach(el => el.classList.remove("open"));
      if (!wasOpen) detailTr.classList.add("open");
    }});

    frag.appendChild(tr);
    frag.appendChild(detailTr);
  }});
  tbody.appendChild(frag);
}}

buildLegend();
buildHead();
buildSourceChips();
["search"].forEach(id => document.getElementById(id).addEventListener("input", scheduleRender));
["qmax", "fcmin"].forEach(id => document.getElementById(id).addEventListener("input", scheduleRender));
["bioactOnly", "secretedOnly"].forEach(id => document.getElementById(id).addEventListener("change", render));
render();
</script>
</body>
</html>
"""


def write_index(comparisons: list[str], species: list[str]) -> Path:
    links = []
    for c in comparisons:
        dirname = c
        links.append(
            f'<li><a href="{dirname}/compound_summary.html">{c}</a> '
            f'<span style="color:var(--ink-muted)">({dirname}/compound_summary.tsv)</span></li>'
        )
    rollup_links = [
        f'<li><a href="all_significant_features_summary_{sp}.html">{sp} rollup</a> '
        f'<span class="muted">(all_significant_features_summary.tsv &mdash; all species)</span></li>'
        for sp in species
    ]
    index_html = f"""<!DOCTYPE html>
<html lang="en">
<head><meta charset="utf-8"><title>Primary differential feature tables</title>
<style>
  body {{ font: 14px/1.5 -apple-system,"Segoe UI",Helvetica,Arial,sans-serif;
         padding: 24px; max-width: 900px; }}
  a {{ color: #0072B2; text-decoration: none; }}
  a:hover {{ text-decoration: underline; }}
  h1 {{ font-size: 18px; }}
  .muted {{ color: #6b6b74; font-size: 13px; }}
  h2 {{ font-size: 14px; margin: 20px 0 6px; }}
</style></head>
<body>
<h1>Primary differential feature tables (8 contrasts, q &lt; 0.05)</h1>
<p class="muted">Each <code>compound_summary.html</code> is a self-contained sortable/filterable
table of significant features with SIRIUS identity. Open directly; no server needed.</p>
<h2>All-significant rollups (chunked per species to stay under GitHub's file-size limit)</h2>
<ul>
  {chr(10).join(rollup_links)}
</ul>
<h2>Per-contrast tables</h2>
<ul>
  {chr(10).join(links)}
</ul>
</body>
</html>
"""
    p = OUT_ROOT / "feature_tables_index.html"
    p.write_text(index_html)
    return p


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--significant", type=Path, default=SIG_TSV)
    args = ap.parse_args()

    sig = Path(args.significant)
    if not sig.exists():
        sys.exit(f"not found: {sig} (run differential_features_primary.py first)")
    df = pd.read_csv(sig, sep="\t")
    if df.empty:
        sys.exit(f"no rows in {sig}")

    comparisons = sorted(df["comparison"].unique())
    written = []
    for comp in comparisons:
        sub = df[df["comparison"] == comp].copy()
        # drop the comparison column for the per-comparison table (it's constant)
        sub_t = sub.drop(columns=["comparison"])
        comp_dir = OUT_ROOT / comp
        comp_dir.mkdir(parents=True, exist_ok=True)

        ordered = [c for c in sub_t.columns]
        sub_t.to_csv(comp_dir / "compound_summary.tsv", sep="\t", index=False)

        title = f"{comp} — significant features"
        html_out = comp_dir / "compound_summary.html"
        html_out.write_text(build_html(sub_t, title))
        written.append((comp, len(sub)))
        print(f"{comp}: {len(sub)} significant rows -> {comp_dir}/compound_summary.html", file=sys.stderr)

    df.to_csv(OUT_ROOT / "all_significant_features_summary.tsv", sep="\t", index=False)
    species_vals = sorted(df["species"].dropna().unique()) if "species" in df.columns else []
    for sp in species_vals:
        spdf = df[df["species"] == sp].copy()
        rollup_html = OUT_ROOT / f"all_significant_features_summary_{sp}.html"
        rollup_html.write_text(
            build_html(spdf, f"All significant features — {sp} (across the 8 primary contrasts)", is_rollup=True)
        )
        print(f"rollup({sp}): {len(spdf)} significant rows -> {rollup_html}", file=sys.stderr)
        del spdf

    idx = write_index(comparisons, species_vals)
    print(f"index: {idx}", file=sys.stderr)


if __name__ == "__main__":
    main()
