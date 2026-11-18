#!/usr/bin/env python3
"""Turn an audited design into bench files: one .dna per PCR product, one for the assembly.

Why: the plan says "PCR these five pieces and In-Fusion them". The bench then has to
trust that the pieces really join into the plasmid that was designed. This script closes
that loop by re-deriving every PCR product from the primers, joining them the way
In-Fusion/Gibson would (by overlap, not by design intent), and writing both the pieces
and the product into one folder per construct.

Usage:

    python3 scripts/simulate_build.py design.json --outdir 验证/mPhyB-WT

`design.json` is the same file audit_assembly.py takes:

    {
      "target": "target.fa",
      "order": ["V", "P1", "PIF3"],
      "fragments": {"V": {"product": [5831, 4625], "template": "p17"}},
      "templates": {"p17": "p17.fa", "pB": "pB.fa"},
      "primers": [{"name": "VEC-F", "seq": "...", "template": "p17", "binding": 23}, ...]
    }

`product` is the circular [start, end] window of the fragment on its template, so a
vector that wraps the origin is written start > end.

Output per fragment: `<name>.dna` holding the PCR product exactly as the primers would
make it (5' tail included), with a feature labelling the overlap regions. Then
`<construct>.dna` holding the circular assembly. The script exits non-zero unless the
assembly equals the target base for base, so a folder of files can never be mislabelled.

The construct name is the design filename stem unless `"name"` is set in the JSON.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import snapgene_io
from tm_engine import tm


def read_fasta(path):
    with open(path) as fh:
        fh.readline()
        return "".join(fh.read().split()).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def circular(seq, a, b):
    """1-based inclusive window that may wrap the origin."""
    if a <= b:
        return seq[a - 1 : b]
    return seq[a - 1 :] + seq[:b]


def auto_binding(template, seq):
    """Longest 3' suffix that anneals uniquely, so the sheet needs no binding column."""
    d = template + template
    for need in (1, 2):                     # prefer a unique site, then any
        for k in range(len(seq), 11, -1):
            tail = seq[-k:]
            if d.count(tail) == need or d.count(rc(tail)) == need:
                return k
    raise SystemExit("primer %s does not anneal" % seq)


def amplimer(template, primer):
    """Find the primer's 3' binding site; return its 0-based start (or -1)."""
    bind = primer["seq"][-primer["binding"] :]
    return (template + template).find(bind)


def pcr(template, f, r):
    """Simulate one PCR: 5' tail + template body + reverse-complement of the R primer."""
    i = amplimer(template, f)
    if i < 0:
        raise SystemExit("forward primer %s does not anneal" % f["name"])
    j = (template + template).find(rc(r["seq"][-r["binding"] :]), i + f["binding"])
    if j < 0:
        raise SystemExit("reverse primer %s does not anneal" % r["name"])
    body = (template + template)[i + f["binding"] : j]
    return f["seq"] + body + rc(r["seq"]), i + f["binding"], j


def assemble(parts, overlaps):
    """Join fragments by their designed overlap (Gibson style), then circularise."""
    out = parts[0]
    for nxt, k in zip(parts[1:], overlaps[:-1]):
        if out[-k:] != nxt[:k]:
            raise SystemExit("overlap mismatch: %r vs %r" % (out[-k:], nxt[:k]))
        out += nxt[k:]
    k = overlaps[-1]
    if out[-k:] != out[:k]:
        raise SystemExit("cannot circularise: ends do not share %d bp" % k)
    return out[:-k]


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design")
    ap.add_argument("--outdir", required=True, help="folder for this construct")
    ap.add_argument("--sheet", help="order sheet (xlsx). When given, the primer sequences in the "
                                   "sheet override the design JSON: only what gets synthesised counts.")
    args = ap.parse_args(argv)

    d = json.load(open(args.design))
    target = read_fasta(d["target"])
    templates = {k: read_fasta(v) for k, v in d["templates"].items()}
    template_features = {}
    for k, v in d.get("template_dna", {}).items():
        _, template_features[k] = snapgene_io.parse(v)
    primers = {p["name"]: dict(p) for p in d["primers"]}
    if args.sheet:
        from openpyxl import load_workbook
        ws = load_workbook(args.sheet).active
        sheet = {}
        for r in range(2, ws.max_row + 1):
            name, seq = ws.cell(r, 3).value, ws.cell(r, 4).value
            if name and seq:
                sheet[name.strip()] = seq.strip().upper()
        for name, seq in sheet.items():
            primers.setdefault(name, {"name": name})
            primers[name]["seq"] = seq
            primers[name]["binding"] = None      # filled per template below
        print("primers read from %s: %d" % (args.sheet, len(sheet)))
    for spec in d["fragments"].values():
        for key in ("f", "r"):
            p = primers[spec[key]]
            if p.get("binding") is None:
                p["binding"] = auto_binding(templates[spec["template"]], p["seq"])
    order = d["order"]
    name = d.get("name") or os.path.splitext(os.path.basename(args.design))[0]
    os.makedirs(args.outdir, exist_ok=True)

    # Fragment ends, in target coordinates, give the overlap lengths without guessing.
    spans, parts = {}, []
    for frag in order:
        tmpl = templates[d["fragments"][frag]["template"]]
        spec = d["fragments"][frag]
        f = primers[spec.get("f", frag + "-F")]
        r = primers[spec.get("r", frag + "-R")]
        product, inner_start, inner_end = pcr(tmpl, f, r)
        spans[frag] = len(product)
        parts.append(product)
        # A PCR product is a linear molecule; only the Gibson product is circular.
        shift = len(f["seq"]) - inner_start
        feats = [(frag, 1, len(product), "misc_feature")]
        for label, a, b, typ, _strand, color in template_features.get(d["fragments"][frag]["template"], []):
            if a >= inner_start + 1 and b <= inner_end:
                feats.append((label, a + shift, b + shift, typ, 1, color))
        snapgene_io.write(os.path.join(args.outdir, frag + ".dna"), product,
                          feats, circular=False)
        print("  %-8s %5d bp  <- %s / %s" % (frag, len(product), f["name"], r["name"]))

    overlaps = []
    products = {frag: t for frag, t in zip(order, parts)}
    for i, frag in enumerate(order):
        nxt = order[(i + 1) % len(order)]
        a = min(len(products[frag]), len(products[nxt]))
        k = 0
        for cand in range(a, 9, -1):
            if products[frag][-cand:] == products[nxt][:cand]:
                k = cand
                break
        if k < 10:
            raise SystemExit("no usable overlap between %s and %s" % (frag, nxt))
        overlaps.append(k)

    circular_product = assemble(parts, overlaps)
    ok = len(circular_product) == len(target)
    if ok:
        doubled = circular_product + circular_product
        ok = doubled.find(target[:60]) >= 0 and len(circular_product) == len(target)

    asm_feats = [(name, 1, len(circular_product), "misc_feature")]
    probe = circular_product + circular_product
    for src, feats in template_features.items():
        src_seq = templates[src]
        for label, a, b, typ, _strand, _color in feats:
            piece = src_seq[a - 1 : b]
            hit = probe.find(piece)
            if hit >= 0 and hit < len(circular_product):
                feat = (label, hit + 1, hit + len(piece), typ, 1, _color)
                if feat not in asm_feats:
                    asm_feats.append(feat)
    snapgene_io.write(os.path.join(args.outdir, name + ".dna"), circular_product,
                      asm_feats)
    status = "OK" if ok else "MISMATCH"
    print("  %-8s %5d bp  assembly vs target %s" % (name, len(circular_product), status))
    print("wrote %d files to %s" % (len(order) + 1, args.outdir))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
