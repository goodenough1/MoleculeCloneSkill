#!/usr/bin/env python3
"""Read and write SnapGene .dna files (sequence + a few labelled features).

Why this module exists: the validation step has to hand the bench openable files, one
per PCR fragment plus the assembled plasmid. A FASTA works for a script, a .dna works
for a human who wants to check the junctions in SnapGene.

File layout (same as SnapGene writes): a sequence of segments, each a 1-byte type, a
4-byte big-endian payload length, then the payload. Type 9 is the header carrying the
"SnapGene" magic, type 0 the raw sequence, type 10 the features XML. Types 3 and 11
(enzyme list, edit history) are optional and are not written here.
"""

from __future__ import annotations

import struct
from xml.sax.saxutils import escape

# 14-byte header payload: "SnapGene" + type_of_sequence (0x0001 = DNA)
# + export_version (15) + import_version (20). Topology is NOT here.
HEADER = b"SnapGene\x00\x01\x00\x0f\x00\x14"
# Block 0 starts with a property-flags byte: bit0 circular, bit1 double-stranded.
FLAGS_CIRCULAR = 0x03
FLAGS_LINEAR = 0x02


def parse(path):
    """Return (sequence, [(label, start, end, type, strand), ...])."""
    with open(path, "rb") as fh:
        data = fh.read()
    pos, seq, xml = 0, None, None
    while pos + 5 <= len(data):
        typ = data[pos]
        ln = struct.unpack(">I", data[pos + 1 : pos + 5])[0]
        payload = data[pos + 5 : pos + 5 + ln]
        if typ == 0:
            # payload[0] is the property-flags byte, not a base
            seq = payload[1:].decode("ascii", "replace")
        elif typ == 10:
            xml = payload.decode("utf-8", "replace")
        pos += 5 + ln
    return (seq or "").upper(), features_from_xml(xml or "")


def features_from_xml(xml):
    import re

    out = []
    for blk in re.findall(r"<Feature\b.*?</Feature>", xml, re.S):
        name = re.search(r'name="([^"]*)"', blk)
        typ = re.search(r'type="([^"]*)"', blk)
        for a, b, color in re.findall(r'<Segment range="(\d+)-(\d+)"(?: color="([^"]*)")?', blk):
            out.append((name.group(1) if name else "", int(a), int(b),
                        typ.group(1) if typ else "misc_feature", 1, color or ""))
    return out


def features_to_xml(features):
    parts = ['<?xml version="1.0"?><Features nextValidID="%d">' % (len(features) + 1)]
    for i, feat in enumerate(features):
        label, start, end, typ = feat[:4]
        strand = feat[4] if len(feat) > 4 else 1
        color = feat[5] if len(feat) > 5 and feat[5] else ("#7f9dc4" if typ == "CDS" else "#a6acb3")
        parts.append(
            '<Feature recentID="%d" name="%s" type="%s" allowSegmentOverlaps="0" '
            'consecutiveTranslationNumbering="1" directionality="%d">'
            '<Segment range="%d-%d" color="%s" type="standard"/>%s</Feature>'
            % (i, escape(label, {'"': "&quot;"}), typ, 1 if strand >= 0 else -1,
               start, end, color, "" if strand >= 0 else "")
        )
    parts.append("</Features>")
    return "".join(parts)


def write(path, seq, features=None, circular=True):
    """Write a .dna openable by SnapGene. `features` = [(label, start, end, type)]."""
    seq = seq.upper()
    flags = FLAGS_CIRCULAR if circular else FLAGS_LINEAR
    segs = [(9, HEADER), (0, bytes([flags]) + seq.encode("ascii"))]
    if features:
        segs.append((10, features_to_xml(features).encode("utf-8")))
    with open(path, "wb") as fh:
        for typ, payload in segs:
            fh.write(bytes([typ]))
            fh.write(struct.pack(">I", len(payload)))
            fh.write(payload)
    return path


def write_fasta(path, seq, name):
    with open(path, "w") as fh:
        fh.write(">%s\n" % name)
        fh.write("\n".join(seq[i : i + 70] for i in range(0, len(seq), 70)) + "\n")
    return path
