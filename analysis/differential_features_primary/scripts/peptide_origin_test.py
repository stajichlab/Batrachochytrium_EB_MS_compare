#!/usr/bin/env python3
"""RETRACTED TEST (2026-09-02) -- kept for the audit trail, do not cite its conclusion.

This script tests whether blank-clearing liq peptides are medium-protein digest
products via proline composition parsed from SIRIUS structure names. **The test
fails its own controls and its conclusion is withdrawn.** Run it only to
reproduce the negative result.

Why it fails:
  * The whole SIRIUS annotation table, parsed identically, is 19.5% Pro+Hyp
    (872 peptides, 5,048 residues). The shortlist values -- Bd 24.3%, Bsal
    16.9% -- are indistinguishable from that baseline (binomial p=0.23, 0.42).
    The ~20% level is a property of which peptides the structure database can
    name, not of these samples.
  * Against each species' own non-blank-clearing (medium) peptides (Bd 19.7%,
    Bsal 18.2%) the shortlist differs at Fisher p=0.27 / 0.75. The statistic
    cannot separate a blank-clearing peptide from a medium peptide.
  * The reference substrate was wrong: tryptone is a digest of WHOLE casein
    (~11.3% Pro), not beta-casein (16.7%). Bd's 24.3% rejects whole casein.
  * The unit of independence is the molecule, not the residue: 19 + 31 parsed
    rows are ~15 and ~24 distinct molecules (some are isotope/duplicate-m/z
    copies), so the binomial overstates precision by roughly the peptide length.

Lesson worth keeping: a compositional statistic computed over database-assigned
identities inherits the database's composition. Its control must be a
same-database baseline, not an external reference proteome.

The surviving, much narrower result uses fragment ions with a built-in negative
control: the hydroxyproline immonium (86.0600) is in 26.3% of Bsal liq peptide
MS2 spectra vs 9.2% of Bd's (OR 3.53, p=7.1e-16). Hyp is collagen-specific and
only Bsal's medium contains gelatin. See CORRECTED_REANALYSIS_REPORT.md 6.2.

---- original docstring below ----

Test whether the blank-clearing 'secreted' peptides are medium-protein digest products.

The question
------------
After the 2026-09-02 corrections (media blanks removed from the analysis
matrix, artifact rows dropped, paired blank filter, exact permutation null,
MS2 gate), the defensible liq-enriched shortlist is small and almost
entirely one compound class: SIRIUS calls 56/90 (Bd) and 102/139 (Bsal)
shortlist-ready features "Amino acids and Peptides", and the
highest-confidence structures are short PROLINE-RICH peptides --
Val-Leu-Pro-Val-Pro, Pro-Val-Val-Pro, Tyr-Pro-Phe-Pro, H-Val-Val-Pro-Pro-Phe.

Two hypotheses explain a liq-enriched, blank-clearing peptide:

  H1 (biosynthesis)  a non-ribosomal peptide made by the tier-1 NRPS.
  H2 (proteolysis)   a fragment of MEDIUM PROTEIN released by secreted
                     fungal proteases. Bd's medium is 1% tryptone (casein
                     digest); Bsal's is 50% TGHL (tryptone/gelatin
                     hydrolysate/lactose). Both substrates are unusually
                     proline-rich: bovine beta-casein is ~17% Pro and
                     collagen/gelatin ~22% (Pro+Hyp), against ~5% Pro in an
                     average proteome.

H2 also has independent genomic support in this very repo: the strongest
comparative-genomics result is a large MEROPS M36 fungalysin expansion in
Bsal (233/247 secreted protease candidates; PROTEASE_CANDIDATES.md), i.e.
the genome predicts exactly this activity.

The test
--------
Parse amino-acid composition out of the SIRIUS structure names of the
blank-clearing, MS2-backed peptide features (only names that are
unambiguous residue strings; everything else is reported as unparsed rather
than guessed), then compare observed Pro frequency against three reference
compositions with a binomial test:

    beta-casein (bovine, P02666)        ~17.0% Pro
    collagen/gelatin (Pro+Hyp)          ~22.0%
    average proteome                     ~5.0% Pro

H2 predicts Pro frequency near the casein/gelatin values and far above the
proteome baseline. H1 predicts no particular Pro enrichment (NRPS A-domain
specificity for the Bd cluster was predicted as Ala, not Pro).

This does not prove H2 for any individual feature -- a definitive answer
needs a defined/labelled medium (see the report's experimental
recommendation). It tests whether the CLASS is dominated by digest products.

Usage: python3 scripts/peptide_origin_test.py
Outputs: lifestage_trend/../peptide_origin/{peptide_origin.tsv, peptide_composition.png}
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest

REPO = Path(__file__).resolve().parents[3]
CUR = REPO / "analysis" / "differential_features_primary" / "liq_enriched_curation"
OUT = REPO / "analysis" / "differential_features_primary" / "peptide_origin"

AA3 = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
       "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val",
       "Hyp", "Nle", "Orn"]

# Reference Pro fractions. Casein/collagen values are standard published
# compositions; the proteome baseline is the canonical ~5% figure.
REFERENCES = {
    "beta-casein (bovine)": 0.170,
    "collagen/gelatin (Pro+Hyp)": 0.220,
    "average proteome": 0.050,
}

# A name is treated as a peptide only if it is (almost) entirely residue
# tokens plus standard peptide decorations -- never inferred from a trade
# name, and never from a formula.
_TOKEN = re.compile("|".join(AA3))
_DECOR = re.compile(r"(^H-|-OH$|-NH2$|^Ac-|^Boc-|^cyclo\[|\]$|[-\(\)\[\],\s]|"
                    r"\bDL\b|\bD\b|\bL\b|xi|\d+)", re.IGNORECASE)


def parse_residues(name: str) -> list[str] | None:
    """Residue list if `name` is an unambiguous residue string, else None."""
    if not isinstance(name, str) or not name.strip():
        return None
    toks = _TOKEN.findall(name)
    if len(toks) < 2:
        return None
    # Reject names where meaningful non-residue text remains: that means the
    # name is a systematic/trivial name that merely happens to contain a
    # 3-letter substring (e.g. "Prodigiosin" contains "Pro").
    residual = _DECOR.sub("", _TOKEN.sub("", name))
    if len(residual) > 2:
        return None
    return toks


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows, per_species = [], {}

    for sp in ["dendrobatidis", "salamandrivorans"]:
        f = CUR / f"{sp}_liq_enriched_top.tsv"
        if not f.exists():
            sys.exit(f"missing {f} -- run `pixi run usi-curation` first")
        d = pd.read_csv(f, sep="\t")
        for col in ("shortlist_ready", "sirius_structure_name"):
            if col not in d.columns:
                sys.exit(f"{f} lacks {col} -- re-run usi-curation")
        ready = d[d["shortlist_ready"] & d["sirius_structure_name"].notna()]

        parsed, unparsed = [], 0
        for _, r in ready.iterrows():
            res = parse_residues(r["sirius_structure_name"])
            if res is None:
                unparsed += 1
                continue
            parsed.append({
                "species": sp, "row_id": int(r["row_id"]), "mz": r["mz"],
                "confidence": r.get("sirius_structure_confidence"),
                "name": r["sirius_structure_name"],
                "n_residues": len(res),
                "n_pro": sum(1 for x in res if x in ("Pro", "Hyp")),
                "residues": ";".join(res),
            })
        pdf = pd.DataFrame(parsed)
        per_species[sp] = (pdf, len(ready), unparsed)
        rows.append(pdf)

        if pdf.empty:
            print(f"{sp}: no parseable peptide names among {len(ready)} ready+named",
                  file=sys.stderr)
            continue
        tot_res = int(pdf["n_residues"].sum())
        tot_pro = int(pdf["n_pro"].sum())
        frac = tot_pro / tot_res
        print(f"\n=== {sp} ===", file=sys.stderr)
        print(f"  shortlist-ready & named: {len(ready)}; parsed as peptides: {len(pdf)} "
              f"(unparsed/non-peptide names: {unparsed})", file=sys.stderr)
        print(f"  residues={tot_res}  Pro(+Hyp)={tot_pro}  frequency={frac:.1%}",
              file=sys.stderr)
        for label, p0 in REFERENCES.items():
            bt = binomtest(tot_pro, tot_res, p0, alternative="two-sided")
            verdict = "consistent" if bt.pvalue > 0.01 else "rejected"
            print(f"    vs {label:28s} (Pro={p0:.1%}): p={bt.pvalue:.2e}  {verdict}",
                  file=sys.stderr)

    allp = pd.concat([r for r in rows if not r.empty], ignore_index=True)
    allp.to_csv(OUT / "peptide_origin.tsv", sep="\t", index=False)

    # Figure: observed Pro frequency vs reference compositions.
    fig, ax = plt.subplots(figsize=(8, 4.6))
    labels, vals, cols, errs = [], [], [], []
    for sp, colour in [("dendrobatidis", "#D55E00"), ("salamandrivorans", "#0072B2")]:
        pdf = per_species[sp][0]
        if pdf.empty:
            continue
        n, k = int(pdf["n_residues"].sum()), int(pdf["n_pro"].sum())
        f = k / n
        labels.append(f"{sp}\n(n={k}/{n} residues)")
        vals.append(f * 100)
        cols.append(colour)
        errs.append(1.96 * np.sqrt(f * (1 - f) / n) * 100)
    ax.bar(range(len(vals)), vals, yerr=errs, capsize=4, color=cols,
           edgecolor="white", width=0.55)
    for i, (v, e) in enumerate(zip(vals, errs)):
        ax.text(i, v + e + 0.7, f"{v:.1f}%", ha="center", fontsize=10, fontweight="bold")
    styles = {"beta-casein (bovine)": ("--", "#444444"),
              "collagen/gelatin (Pro+Hyp)": (":", "#777777"),
              "average proteome": ("-.", "#B00020")}
    for label, p0 in REFERENCES.items():
        ls, c = styles[label]
        ax.axhline(p0 * 100, linestyle=ls, color=c, linewidth=1.3,
                   label=f"{label}  {p0:.0%}")
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylabel("Pro (+Hyp) frequency  (%)")
    ax.set_title("Blank-clearing, MS2-backed liq peptides are proline-rich,\n"
                 "matching casein/gelatin digest rather than an average proteome",
                 fontsize=10)
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "peptide_composition.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "peptide_composition.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"\nwrote {len(allp)} parsed peptides to {OUT}", file=sys.stderr)


if __name__ == "__main__":
    main()
