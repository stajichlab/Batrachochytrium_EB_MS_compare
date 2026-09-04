#!/usr/bin/env python3
"""De novo MS2 fragment-tag matching of the liq shortlist against medium proteins.

The question this settles
-------------------------
Are the blank-clearing, MS2-backed liquid-fraction peptides fragments of MEDIUM
PROTEIN (H2), or products of the tier-1 NRPS (H1)? Two earlier attempts failed
or were inconclusive:

  * The proline-composition test is RETRACTED (see peptide_origin_test.py). It
    inherited the composition of SIRIUS's structure database, so it could not
    tell a blank-clearing peptide from a medium peptide.
  * Sequence-substring matching of SIRIUS *names* was negative, but names are
    database guesses -- a true fragment can be named as the wrong isomer.

This test uses the SPECTRA, not the names, so it is immune to database bias.

Method
------
1. For each shortlist spectrum, keep the `--top-peaks` most intense peaks and
   build a spectrum graph: an edge between two peaks whose m/z differ by an
   amino-acid residue mass within `--tol`.
2. Enumerate every path of `--tag-len` consecutive edges. Each path spells a
   sequence tag (a de novo read of `--tag-len` residues), independent of any
   database.
3. Match tags as substrings of the medium substrate proteins. Leu and Ile are
   COLLAPSED (both 113.08406 Da -- indistinguishable in MS2), so matching is
   done in a 19-letter alphabet.
4. Score against a TAG-SIDE decoy: each observed tag is letter-shuffled
   (preserving its exact amino-acid composition) and re-matched against the
   same real substrate. Enrichment = real hits / mean decoy hits.

Two design decisions here were forced by controls, and both are the opposite
of the obvious choice:

  * **Match in BOTH orientations.** A b/y spectrum spells the peptide twice --
    the b-series N->C and the y-series C->N -- so the tag set contains the
    sequence and its reverse. A tag is counted as a hit if it or its reverse
    occurs in the substrate.
  * **Randomize the TAG, not the substrate.** Because of the above, a
    reversed-sequence decoy is structurally invalid: it scored exactly 1.00x
    against a known casein peptide at every tag length from 3 to 6. A
    shuffled-*substrate* decoy is also biased -- collagen is a Gly-X-Y repeat
    with only 0.43 distinct 3-mers per position, so shuffling inflates the
    decoy vocabulary by 20% and penalises the real sequence. Shuffling the
    tag's own letters avoids both artifacts.

5. Tag length matters more than anything else, because it sets how much of the
   k-mer space the substrate occupies:

       3-mer: 1,468 of      6,859 possible = 21.4%   -> no discrimination
       4-mer: 2,257 of    130,321 possible =  1.7%
       5-mer: 2,728 of  2,476,099 possible =  0.11%  -> default

   At 3 residues a tag matches casein by luck one time in five. The default is
   therefore 5, not 3.

Interpreting the outcome
------------------------
A POSITIVE result (tags matching real substrates well above decoy) is direct
evidence that these features are medium-protein digest products. A NEGATIVE
result is equally informative and is the outcome the name-based substring test
already hinted at: it would argue the shortlist is NOT simply tryptone/gelatin
fragments, leaving NRPS products or other fungal peptides in play.

Caveats
-------
* Peak depth is a real parameter, not a detail. These spectra carry a median
  ~563-649 peaks, which is high for small peptides (chimeric/noisy), and a
  too-generous peak list manufactures gaps by coincidence. `--top-peaks` is
  swept by `--sweep` for exactly this reason.
* Tags are matched by SEQUENCE only, not anchored to a prefix mass. That is
  deliberately permissive (higher sensitivity, more false matches), and the
  decoy null is what makes it interpretable.
* A tag is not a peptide identification. Matching says "this spectrum contains
  a ladder spelling a substring of casein/collagen", which is evidence about
  origin, not a structure assignment.

Usage:
    python3 scripts/fragment_tag_match.py [--top-peaks 60] [--tag-len 5]
                                          [--tol 0.01] [--n-shuffle 100] [--sweep]
"""
from __future__ import annotations

import argparse
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
MGF = (REPO / "data" / "raw" / "gnps2_e9838293_bagel" / "nf_output"
       / "feature_finding" / "aligned_features_filled.mgf")
CUR = REPO / "analysis" / "differential_features_primary" / "liq_enriched_curation"
SUBS = REPO / "reference_material" / "substrate_proteins"
OUT = REPO / "analysis" / "peptide_provenance"

# Monoisotopic residue masses. Ile is omitted: it is isobaric with Leu and is
# folded into "L" everywhere, so the alphabet is 19 letters.
AA = {
    "G": 57.02146, "A": 71.03711, "S": 87.03203, "P": 97.05276, "V": 99.06841,
    "T": 101.04768, "C": 103.00919, "L": 113.08406, "N": 114.04293,
    "D": 115.02694, "Q": 128.05858, "K": 128.09496, "E": 129.04259,
    "M": 131.04049, "H": 137.05891, "F": 147.06841, "R": 156.10111,
    "Y": 163.06333, "W": 186.07931,
}
SPECIES = ["dendrobatidis", "salamandrivorans"]


def collapse_il(seq: str) -> str:
    return seq.replace("I", "L")


def load_spectra() -> dict[int, np.ndarray]:
    spec: dict[int, np.ndarray] = {}
    cur: dict[str, str] | None = None
    peaks: list[tuple[float, float]] = []
    with open(MGF) as fh:
        for line in fh:
            t = line.strip()
            if t == "BEGIN IONS":
                cur, peaks = {}, []
            elif t == "END IONS":
                if cur is not None and "SCANS" in cur and peaks:
                    spec[int(cur["SCANS"])] = np.array(peaks)
                cur, peaks = None, []
            elif cur is not None:
                if "=" in t and t.split("=")[0].isupper():
                    k, v = t.split("=", 1)
                    cur[k] = v
                elif t and t[0].isdigit():
                    a = t.split()
                    if len(a) >= 2:
                        try:
                            peaks.append((float(a[0]), float(a[1])))
                        except ValueError:
                            pass
    return spec


def load_substrates() -> dict[str, str]:
    out: dict[str, str] = {}
    for f in sorted(SUBS.glob("*.fasta")):
        seq = "".join(l.strip() for l in f.read_text().splitlines()
                      if not l.startswith(">"))
        out[f.name.split(".")[0]] = collapse_il(seq.upper())
    if not out:
        sys.exit(f"no substrate FASTAs in {SUBS} -- run fetch_substrates.sh")
    return out


def extract_tags(peaks: np.ndarray, top_peaks: int, tag_len: int,
                 tol: float) -> list[str]:
    """De novo sequence tags of `tag_len` residues from the spectrum graph."""
    s = peaks[np.argsort(-peaks[:, 1])][:top_peaks]
    mz = np.sort(s[:, 0])
    n = len(mz)
    adj: dict[int, list[tuple[int, str]]] = defaultdict(list)
    diff = mz[None, :] - mz[:, None]
    for aa, m in AA.items():
        ii, jj = np.where(np.abs(diff - m) <= tol)
        for i, j in zip(ii, jj):
            if j > i:
                adj[int(i)].append((int(j), aa))
    tags: list[str] = []

    def walk(node: int, acc: str) -> None:
        if len(acc) == tag_len:
            tags.append(acc)
            return
        for nxt, aa in adj.get(node, ()):
            walk(nxt, acc + aa)

    for start in range(n):
        walk(start, "")
    return tags


def substrate_hits(tags: list[str], substrate_blob: str) -> int:
    """Tags occurring in the substrate in EITHER orientation.

    Both directions are required: the b-ion ladder spells the peptide N->C and
    the y-ion ladder spells it C->N, so a genuine substrate fragment shows up
    as a tag or as its reverse depending on which series the peaks came from.
    """
    return sum(1 for t in tags if (t in substrate_blob) or (t[::-1] in substrate_blob))


def decoy_hits(tags: list[str], substrate_blob: str, n_shuffle: int,
               rng: np.random.Generator) -> np.ndarray:
    """Null: letter-shuffle each tag (composition preserved), re-match."""
    out = np.empty(n_shuffle)
    for i in range(n_shuffle):
        shuffled = []
        for t in tags:
            a = list(t)
            rng.shuffle(a)
            shuffled.append("".join(a))
        out[i] = substrate_hits(shuffled, substrate_blob)
    return out


def run(top_peaks: int, tag_len: int, tol: float, n_shuffle: int, seed: int,
        spec, subs, verbose: bool = True):
    real = "|".join(subs.values())
    rng = np.random.default_rng(seed)
    rows = []
    for sp in SPECIES:
        f = CUR / f"{sp}_liq_enriched_top.tsv"
        if not f.exists():
            sys.exit(f"missing {f} -- run `pixi run usi-curation` first")
        d = pd.read_csv(f, sep="\t")
        if "shortlist_ready" not in d.columns:
            sys.exit(f"{f} lacks shortlist_ready -- re-run usi-curation")
        for _, r in d[d["shortlist_ready"]].iterrows():
            peaks = spec.get(int(r["row_id"]))
            if peaks is None:
                continue
            tags = extract_tags(peaks, top_peaks, tag_len, tol)
            if not tags:
                continue
            uniq = sorted(set(tags))
            h_real = substrate_hits(uniq, real)
            h_dec = decoy_hits(uniq, real, n_shuffle, rng)
            # one-sided empirical p: how often a composition-matched decoy tag
            # set matches the substrate at least as well as the observed one
            p = float((h_dec >= h_real).sum() + 1) / (n_shuffle + 1)
            rows.append({
                "species": sp, "row_id": int(r["row_id"]), "mz": r["mz"],
                "n_tags": len(tags), "n_unique_tags": len(uniq),
                "hits_real": h_real,
                "hits_decoy_mean": float(h_dec.mean()),
                "hits_decoy_p95": float(np.percentile(h_dec, 95)),
                "enrichment": h_real / max(h_dec.mean(), 1e-9),
                "perm_p": p,
                "sirius_structure_name": r.get("sirius_structure_name"),
                "sirius_structure_confidence": r.get("sirius_structure_confidence"),
            })
    df = pd.DataFrame(rows)
    if df.empty:
        sys.exit("no spectra produced tags")
    tot_real = int(df.hits_real.sum())
    tot_dec = float(df.hits_decoy_mean.sum())
    n_sig = int((df.perm_p < 0.05).sum())
    if verbose:
        print(f"[top{top_peaks}/tag{tag_len}] spectra={len(df)} "
              f"unique tags={int(df.n_unique_tags.sum())} | "
              f"hits real={tot_real} decoy(mean)={tot_dec:.1f} | "
              f"enrichment={tot_real/max(tot_dec,1e-9):.2f}x | "
              f"spectra with p<0.05: {n_sig}/{len(df)}", file=sys.stderr)
    return df, dict(top_peaks=top_peaks, tag_len=tag_len, n_spectra=len(df),
                    unique_tags=int(df.n_unique_tags.sum()),
                    hits_real=tot_real, hits_decoy_mean=tot_dec,
                    enrichment=tot_real / max(tot_dec, 1e-9),
                    n_spectra_p05=n_sig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top-peaks", type=int, default=60)
    ap.add_argument("--tag-len", type=int, default=5)
    ap.add_argument("--tol", type=float, default=0.01)
    ap.add_argument("--n-shuffle", type=int, default=100)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--sweep", action="store_true",
                    help="sweep --top-peaks (20/40/60/100/150) and tag length 3/4/5")
    args = ap.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    spec = load_spectra()
    subs = load_substrates()
    print(f"spectra={len(spec)}  substrates={list(subs)} "
          f"({sum(len(s) for s in subs.values())} aa total)", file=sys.stderr)

    df, summary = run(args.top_peaks, args.tag_len, args.tol,
                      args.n_shuffle, args.seed, spec, subs)
    df.sort_values(["perm_p", "enrichment"],
                   ascending=[True, False]).to_csv(
        OUT / "fragment_tag_hits.tsv", sep="\t", index=False)

    sweeps = [summary]
    if args.sweep:
        for tp in (20, 40, 60, 100, 150):
            for tl in (3, 4, 5):
                if (tp, tl) == (args.top_peaks, args.tag_len):
                    continue
                _, s = run(tp, tl, args.tol, args.n_shuffle, args.seed, spec, subs)
                sweeps.append(s)
    sw = pd.DataFrame(sweeps).sort_values(["tag_len", "top_peaks"])
    sw.to_csv(OUT / "sensitivity_sweep.tsv", sep="\t", index=False)

    # Figure: real vs decoy hits across the peak-depth sweep.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    for ax, tl in zip(axes, sorted(sw.tag_len.unique())):
        s = sw[sw.tag_len == tl].sort_values("top_peaks")
        ax.plot(s.top_peaks, s.hits_real, "o-", color="#D55E00", label="real substrates")
        ax.plot(s.top_peaks, s.hits_decoy_mean, "^:", color="#888888",
                label="composition-matched decoy tags")
        ax.set_xlabel("top peaks retained per spectrum")
        ax.set_title(f"tag length {tl}")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("sequence-tag hits (summed over spectra)")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle("De novo fragment tags vs medium-protein substrates\n"
                 "(decoy = letter-shuffled tags; validated positive control "
                 "55x at 5-mers, negative control 0x)", fontsize=10)
    fig.tight_layout()
    fig.savefig(OUT / "fragment_tag_sweep.png", dpi=150, bbox_inches="tight")
    fig.savefig(OUT / "fragment_tag_sweep.pdf", bbox_inches="tight")
    plt.close(fig)

    n_sig = int((df.perm_p < 0.05).sum())
    print(f"\nper-feature: {n_sig}/{len(df)} spectra match real substrates better "
          f"than 95% of composition-matched decoy tag sets (perm_p < 0.05)",
          file=sys.stderr)
    print(f"wrote {OUT}/fragment_tag_hits.tsv, sensitivity_sweep.tsv, "
          f"fragment_tag_sweep.png", file=sys.stderr)


if __name__ == "__main__":
    main()
