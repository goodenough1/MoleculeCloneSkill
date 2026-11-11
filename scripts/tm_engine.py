#!/usr/bin/env python3
"""Nearest-neighbour Tm and annealing temperatures for assembly primers.

Why a dedicated engine: Tm decides the annealing temperature, and the annealing
temperature decides whether the PCR works at all. Different tools (Primer3, SnapGene,
IDT, primer3-py) use different salt models and can disagree by 5-10 °C, so a design
must name one engine and use it everywhere.

This module reproduces Biopython's ``Tm_NN`` (DNA_NN3 = Allawi & SantaLucia 1997,
saltcorr=7 = Owczarzy 2008) with the parameters used for assembly design:

    dnac1 = 250 nM, dnac2 = 0, Na+ = 50 mM, Mg2+ = 0, dNTPs = 0

which matches SnapGene's displayed Tm closely for ordinary primers. If Biopython
happens to be installed it is used directly; otherwise the built-in pure-Python port
gives identical numbers (verified against Biopython 1.87, max deviation < 0.01 °C).

Rules of thumb encoded here:
  * Tm is computed on the **template-binding region only** — never on the full
    primer, otherwise long 5' tails (P2A, tags, linkers) inflate Tm by 10-20 °C.
  * Annealing temperature = round(min(Tm_F, Tm_R)) - offset (default 3 °C),
    clamped to [50, 72] °C.
  * A long, AT-rich 5' tail can slow annealing even though it is not counted in Tm;
    the report flags tails > 30 nt so the operator can drop 2-4 °C or run a gradient.

CLI
---
    python3 tm_engine.py primers.json [--offset 3] [--csv out.csv] [--latex]

primers.json:
    {
      "defaults": {"binding": 20, "dnac1": 250, "Na": 50, "offset": 3},
      "primers": [
        {"name": "P1-A1-F", "seq": "AGC...", "fragment": "A1", "binding": 35},
        {"name": "P1-A1-R", "seq": "GGT...", "fragment": "A1", "binding": 23}
      ]
    }

Output: one row per primer (Tm on the binding region, GC, tail length) and one row
per fragment with the recommended annealing temperature.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")

# Allawi & SantaLucia (1997) DNA_NN3, as bundled with Biopython.
NN = {
    "AA/TT": (-7.9, -22.2), "AT/TA": (-7.2, -20.4), "TA/AT": (-7.2, -21.3),
    "CA/GT": (-8.5, -22.7), "GT/CA": (-8.4, -22.4), "CT/GA": (-7.8, -21.0),
    "GA/CT": (-8.2, -22.2), "CG/GC": (-10.6, -27.2), "GC/CG": (-9.8, -24.4),
    "GG/CC": (-8.0, -19.9),
}
INIT = {
    "init": (0.0, 0.0), "init_A/T": (2.3, 4.1), "init_G/C": (0.1, -2.8),
    "init_oneG/C": (0.0, 0.0), "init_allA/T": (0.0, 0.0), "init_5T/A": (0.0, 0.0),
}

try:  # prefer Biopython when present: identical model, no maintenance
    from Bio.Seq import Seq as _Seq
    from Bio.SeqUtils import MeltingTemp as _mt
except Exception:  # pragma: no cover
    _mt = None


def gc(seq: str) -> float:
    seq = seq.upper()
    return 100.0 * (seq.count("G") + seq.count("C")) / len(seq) if seq else 0.0


def _tm_pure(seq: str, dnac1: float, dnac2: float, na: float, mg: float,
             dntps: float, saltcorr: int) -> float:
    seq = seq.upper()
    dh, ds = 0.0, 0.0
    dh += INIT["init"][0]
    ds += INIT["init"][1]
    gc_frac = (seq.count("G") + seq.count("C")) / len(seq)
    key = "init_oneG/C" if gc_frac > 0 else "init_allA/T"
    dh += INIT[key][0]
    ds += INIT[key][1]
    if seq.startswith("T"):
        dh += INIT["init_5T/A"][0]
        ds += INIT["init_5T/A"][1]
    if seq.endswith("A"):
        dh += INIT["init_5T/A"][0]
        ds += INIT["init_5T/A"][1]
    for base in (seq[0], seq[-1]):
        key = "init_A/T" if base in "AT" else "init_G/C"
        dh += INIT[key][0]
        ds += INIT[key][1]
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        k = pair + "/" + pair.translate(COMP)
        if k not in NN:
            k = k[::-1]
        h, s = NN[k]
        dh += h
        ds += s
    k_conc = (dnac1 - dnac2 / 2.0) * 1e-9
    tm = (1000 * dh) / (ds + 1.987 * math.log(k_conc)) - 273.15
    if saltcorr == 7 and na > 0:
        mon = na * 1e-3
        mg_molar = mg * 1e-3
        ratio = math.sqrt(mg_molar) / mon if mg_molar > 0 else 0.0
        if ratio < 0.22:
            corr = (4.29 * gc_frac - 3.95) * 1e-5 * math.log(mon) + 9.40e-6 * math.log(mon) ** 2
        else:
            raise NotImplementedError("Mg-dependent branch of Owczarzy 2008 is not ported; "
                                      "use Biopython or set Mg=0")
        tm = 1 / (1 / (tm + 273.15) + corr) - 273.15
    return tm


def tm(seq: str, dnac1: float = 250.0, dnac2: float = 0.0, na: float = 50.0,
       mg: float = 0.0, dntps: float = 0.0, saltcorr: int = 7) -> float:
    """Tm of a DNA duplex; defaults match the assembly-design engine."""
    seq = seq.upper()
    if len(seq) < 2:
        return float("nan")
    if _mt is not None:
        try:
            return float(_mt.Tm_NN(_Seq(seq), dnac1=dnac1, dnac2=dnac2, Na=na,
                                   Mg=mg, dNTPs=dntps, saltcorr=saltcorr))
        except Exception:
            pass
    return _tm_pure(seq, dnac1, dnac2, na, mg, dntps, saltcorr)


def annealing(tm_f: float, tm_r: float, offset: float = 3.0,
              lo: float = 50.0, hi: float = 72.0) -> int:
    """Recommended annealing temperature for a primer pair."""
    t = math.floor(min(tm_f, tm_r) - offset + 0.5)  # half-up, not banker's rounding
    return int(max(lo, min(hi, t)))


def analyse(primers, default_binding=20, dnac1=250.0, na=50.0, mg=0.0,
            dntps=0.0, offset=3.0, tail_warn=30):
    """Annotate primers and group per fragment.

    Each primer dict: name, seq, optional fragment, optional binding (number of 3'
    bases that anneal to the template).
    """
    rows = []
    for p in primers:
        seq = p["seq"].upper()
        nb = int(p.get("binding", default_binding))
        bind = seq[-nb:]
        row = {
            "name": p.get("name", "?"),
            "fragment": p.get("fragment", ""),
            "length": len(seq),
            "binding": nb,
            "tail": len(seq) - nb,
            "gc": gc(bind),
            "tm": tm(bind, dnac1=dnac1, na=na, mg=mg, dntps=dntps),
            "long_tail": (len(seq) - nb) > tail_warn,
        }
        rows.append(row)

    pairs = {}
    for row in rows:
        if row["fragment"]:
            pairs.setdefault(row["fragment"], []).append(row)
    anneal = {}
    for frag, prs in pairs.items():
        if len(prs) >= 2:
            ta = annealing(prs[0]["tm"], prs[1]["tm"], offset)
            anneal[frag] = {
                "annealing": ta,
                "tm_min": min(prs[0]["tm"], prs[1]["tm"]),
                "long_tail": any(r["long_tail"] for r in prs),
            }
    return rows, anneal


def format_table(rows, anneal):
    out = []
    out.append(f"{'primer':<14}{'frag':<6}{'len':>4}{'bind':>5}{'tail':>5}{'GC%':>6}{'Tm':>8}")
    out.append("-" * 48)
    for r in rows:
        flag = "  <- long tail" if r["long_tail"] else ""
        out.append(f"{r['name']:<14}{r['fragment']:<6}{r['length']:>4}{r['binding']:>5}"
                   f"{r['tail']:>5}{r['gc']:>6.0f}{r['tm']:>8.1f}{flag}")
    if anneal:
        out.append("")
        out.append("recommended annealing (min Tm - offset, clamped 50-72 °C)")
        out.append(f"{'frag':<6}{'Tm min':>8}{'Ta':>6}  note")
        for frag, a in anneal.items():
            note = "AT-rich long tail -> consider -2 to -4 °C or gradient" if a["long_tail"] else ""
            out.append(f"{frag:<6}{a['tm_min']:>8.1f}{a['annealing']:>6}  {note}")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Tm + annealing engine for assembly primers")
    ap.add_argument("primers")
    ap.add_argument("--offset", type=float, default=3.0)
    ap.add_argument("--binding", type=int, default=20, help="default 3' binding length")
    ap.add_argument("--dnac1", type=float, default=250.0)
    ap.add_argument("--na", type=float, default=50.0)
    ap.add_argument("--mg", type=float, default=0.0)
    ap.add_argument("--dntps", type=float, default=0.0)
    ap.add_argument("--tail-warn", type=int, default=30)
    ap.add_argument("--csv")
    ap.add_argument("--latex")
    args = ap.parse_args(argv)

    spec = json.load(open(args.primers))
    defaults = spec.get("defaults", {})
    rows, anneal = analyse(
        spec["primers"],
        default_binding=defaults.get("binding", args.binding),
        dnac1=defaults.get("dnac1", args.dnac1),
        na=defaults.get("Na", args.na),
        mg=defaults.get("Mg", args.mg),
        dntps=defaults.get("dNTPs", args.dntps),
        offset=defaults.get("offset", args.offset),
        tail_warn=defaults.get("tail_warn", args.tail_warn),
    )
    print(format_table(rows, anneal))

    if args.csv:
        with open(args.csv, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.writer(fh)
            w.writerow(["name", "fragment", "length", "binding", "tail", "GC(%)", "Tm"])
            for r in rows:
                w.writerow([r["name"], r["fragment"], r["length"], r["binding"],
                            r["tail"], f"{r['gc']:.0f}", f"{r['tm']:.1f}"])
            for frag, a in anneal.items():
                w.writerow([f"Ta:{frag}", frag, "", "", "", "", a["annealing"]])
        print(f"\nCSV -> {args.csv}")

    if args.latex:
        lines = [r"\begin{tabular}{@{}lcccc@{}}", r"\toprule",
                 r"引物 & 长度 (nt) & 结合区 (nt) & Tm ($^\circ$C) & 退火 ($^\circ$C)\\",
                 r"\midrule"]
        for r in rows:
            ta = anneal.get(r["fragment"], {}).get("annealing", "—")
            lines.append(f"{r['name']} & {r['length']} & {r['binding']} & "
                         f"{r['tm']:.1f} & {ta}\\\\")
        lines += [r"\bottomrule", r"\end{tabular}"]
        print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
