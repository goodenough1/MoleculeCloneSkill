#!/usr/bin/env python3
"""Hard check: every primer must have exactly one priming-competent binding site.

"Unique" is not "the full sequence appears once". A primer whose 3' half matches a
second place with one or two mismatches will still prime there, giving mixed peaks
(Sanger) or extra bands (colony PCR). This script enforces:

  1. the full primer has exactly ONE perfect match on either strand of the template
     (circular templates are searched across the origin);
  2. no other site has <= --max-mismatch mismatches, unless ALL mismatches fall in
     the last --terminal bases (such a site cannot be extended);
  3. for screening primers, `absent_on`: a list of templates on which the primer
     must have NO site at all (the parental plasmid, every PCR fragment, every
     intermediate). `control` + `control_absent` is the one-template shorthand;
  4. optional "pairs": for a screening primer pair, compute the amplicon on the
     target template and on every template in `absent_on` (same shorthand as
     above). A pair passes only when the target gives exactly one product, its
     size equals `expect`, and NO other template gives any product — i.e. only
     the fully assembled plasmid produces a band. Mark the insert-side primer with
     `absent_on`; the shared primer needs no flag;
  5. templates listed in `linear` are treated as linear PCR products, so products
     that would only close across the sequence ends are not counted;
  6. optionally reports Tm of the binding region so annealing can be set.

Exit code is non-zero if any primer or pair fails.

Input JSON
----------
{
  "templates": {"pN1-018": "pN1-018.fa", "S2": "S2.fa", "A1": "A1.fa"},
  "linear": ["A1"],
  "primers": [
    {"name": "F01", "seq": "TCTATGGCTTCTGAGGCGGAAAGA", "template": "pN1-018"},
    {"name": "S01", "seq": "GAGCCTATGGAAAAACGCCAGCA", "template": "S2", "binding": 23},
    {"name": "S2-CHK-F", "seq": "GGGGGAAACGCCTGGTATC", "template": "S2"},
    {"name": "S2-CHK-R", "seq": "CTTCCAG...", "template": "S2",
     "absent_on": ["pN1-16", "A1", "B1", "C1"]}
  ],
  "pairs": [
    {"name": "S2-CHK", "f": "S2-CHK-F", "r": "S2-CHK-R",
     "template": "S2", "absent_on": ["pN1-16", "A1", "B1"], "expect": 3017}
  ]
}

CLI
---
    python3 check_primers.py primers.json [--max-mismatch 2] [--terminal 2] [--csv out.csv]
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tm_engine import tm, gc  # noqa: E402

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
K = 8  # seed length for mismatch screening


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: str) -> str:
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith(">"):
                parts.append(line.upper())
    if not parts:
        raise SystemExit(f"empty FASTA: {path}")
    return "".join(parts)


def build_index(seq: str):
    """8-mer index of the plus strand of a circular sequence."""
    idx = {}
    for i in range(len(seq)):
        idx.setdefault(seq[i:i + K], []).append(i)
    return idx


def candidate_positions(seq, probe, idx):
    n = len(seq)
    out = set()
    for off in (0, len(probe) // 2 - K // 2, len(probe) - K):
        off = max(0, min(off, len(probe) - K))
        for pos in idx.get(probe[off:off + K], []):
            out.add((pos - off) % n)
    return out


def scan(seq, primer, idx, max_mismatch=2, circular=True):
    """Return list of (strand, pos1based, n_mismatch, mismatch_indices).

    circular=True searches across the origin; linear templates (PCR products)
    stop at the ends.
    """
    n = len(seq)
    d = seq + seq if circular else seq
    hits = []
    for strand, probe in (("+", primer), ("-", rc(primer))):
        for pos in candidate_positions(seq, probe, idx):
            if not circular and pos + len(probe) > n:
                continue
            window = d[pos:pos + len(probe)]
            if len(window) < len(probe):
                continue
            diffs = [k for k, (a, b) in enumerate(zip(window, probe)) if a != b]
            if len(diffs) <= max_mismatch:
                hits.append((strand, pos + 1, len(diffs), diffs))
    hits.sort(key=lambda h: (h[2], h[1]))
    return hits


def three_prime_sites(seq, primer, circular=True, min_run=6, tm_limit=42.0):
    """Sites whose 3' end anneals perfectly for long enough to prime above tm_limit.

    A site can carry many mismatches at its 5' end and still be extended: only the
    3' terminal stretch has to anneal. Counting total mismatches misses those, so
    measure the perfect 3' run itself and its Tm instead.
    """
    n = len(seq)
    d = seq + seq if circular else seq
    out = []
    for strand, probe in (("+", primer), ("-", rc(primer))):
        for pos in range(n):
            if not circular and pos + len(probe) > n:
                continue
            window = d[pos:pos + len(probe)]
            if len(window) < len(probe):
                continue
            run = 0
            while run < len(probe) and window[run] == probe[run]:
                run += 1
            if run < min_run:
                continue
            t = tm(primer[-run:])
            if t > tm_limit:
                diffs = [k for k, (a, b) in enumerate(zip(window, probe)) if a != b]
                out.append((strand, pos + 1, len(diffs), run, round(t, 1)))
    out.sort(key=lambda h: -h[4])
    return out


def check(seq, primer, idx, max_mismatch=2, terminal=2, circular=True, tm_limit=42.0):
    hits = scan(seq, primer, idx, max_mismatch, circular)
    perfect = [h for h in hits if h[2] == 0]
    risky = [h for h in hits if h[2] > 0
             and not all(k >= len(primer) - terminal for k in h[3])]
    hot = [h for h in three_prime_sites(seq, primer, circular, tm_limit=tm_limit)
           if h[2] > 0]
    ok = len(perfect) == 1 and not risky and not hot
    return ok, perfect, risky, hits, hot


def perfect_sites(seq, primer, idx, max_mismatch=2, strand=None, circular=True):
    """Positions (1-based, plus strand) of perfect primer matches, optionally by strand."""
    hits = scan(seq, primer, idx, max_mismatch, circular)
    pos = [h[1] for h in hits if h[2] == 0 and (strand is None or h[0] == strand)]
    return sorted(set(pos))


def pair_products(seq, idx, fseq, rseq, max_mismatch=2, circular=True):
    """Amplicon sizes for a primer pair on a circular template.

    F anneals where the plus strand reads fseq (hit on '+'), R where it reads
    rc(rseq) (hit on '-'). A product needs the R site downstream of F in the
    extension direction; anything outside [len(F)+len(R), len(seq)] is not a real
    amplicon.
    """
    n = len(seq)
    fs = perfect_sites(seq, fseq, idx, max_mismatch, strand="+", circular=circular)
    rs = perfect_sites(seq, rseq, idx, max_mismatch, strand="-", circular=circular)
    sizes = set()
    for p in fs:
        for q in rs:
            if circular:
                size = ((q - p) % n) + len(rseq)
            else:
                if q <= p or q + len(rseq) - 1 > n:
                    continue
                size = q - p + len(rseq)
            if len(fseq) + len(rseq) <= size <= n:
                sizes.add(size)
    return sorted(sizes), fs, rs


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("primers")
    ap.add_argument("--max-mismatch", type=int, default=2)
    ap.add_argument("--terminal", type=int, default=2,
                    help="mismatches allowed only within this many 3'-terminal bases")
    ap.add_argument("--binding", type=int, default=None,
                    help="3' bases that anneal; default = whole primer (sequencing primers)")
    ap.add_argument("--csv")
    args = ap.parse_args(argv)

    spec = json.load(open(args.primers))
    tmpl_paths = spec.get("templates", {})
    linear = set(spec.get("linear") or [])
    cache = {}
    failed = 0
    rows = []
    for p in spec["primers"]:
        name = p.get("name", "?")
        seq = p["seq"].upper()
        tname = p.get("template")
        nb = int(p["binding"]) if "binding" in p else (args.binding or len(seq))
        bind = seq[-nb:]
        if tname not in tmpl_paths:
            print(f"[SKIP] {name}: no template '{tname}'")
            continue
        if tname not in cache:
            t = read_fasta(tmpl_paths[tname])
            cache[tname] = (t, build_index(t), tname not in linear)
        t, idx, circ = cache[tname]
        ok, perfect, risky, hits, hot = check(t, bind, idx, args.max_mismatch, args.terminal, circ)
        absent = list(p.get("absent_on") or [])
        if p.get("control") and p.get("control_absent"):
            absent.append(p["control"])
        absent = [x for x in dict.fromkeys(absent) if x != tname]
        clean = []
        control_fail = ""
        for cname in absent:
            if cname not in tmpl_paths:
                ok = False
                control_fail = f"对照 {cname} 没有给出序列文件"
                continue
            if cname not in cache:
                ct = read_fasta(tmpl_paths[cname])
                cache[cname] = (ct, build_index(ct), cname not in linear)
            ct, cidx, ccirc = cache[cname]
            c_ok, c_perfect, c_risky, _, c_hot = check(ct, bind, cidx, args.max_mismatch, args.terminal, ccirc)
            if c_perfect or c_risky or c_hot:
                ok = False
                control_fail = (f"对照 {cname} 有 {len(c_perfect)} 个完全匹配、"
                                f"{len(c_risky)} 个可延伸近似位点（这条引物不是区分侧，去掉 absent_on）")
            else:
                clean.append(cname)
        parts = []
        if clean:
            parts.append(f"{'、'.join(clean)} 无位点")
        if control_fail:
            parts.append(control_fail)
        control_note = ("；" + "；".join(parts)) if parts else ""
        rows.append(dict(name=name, template=tname, length=len(seq),
                         tm=round(tm(bind), 1), gc=round(gc(bind)),
                         exact=len(perfect), risky=len(risky),
                         hot="; ".join(f"{h[0]}{h[1]}({h[3]}nt@{h[4]})" for h in hot),
                         pos=perfect[0][1] if perfect else "",
                         other="; ".join(f"{h[0]}{h[1]}({h[2]}mm@{h[3]})" for h in risky) + control_note))
        if not ok:
            failed += 1
            detail = []
            if len(perfect) != 1:
                detail.append(f"完整匹配 {len(perfect)} 处")
            for h in risky:
                detail.append(f"模板{h[0]}链 {h[1]:,} 位，{h[2]} 个错配（位置 {h[3]}）")
            for h in hot:
                detail.append(f"模板{h[0]}链 {h[1]:,} 位，3′ 端 {h[3]} nt 完美互补、Tm {h[4]} °C，会延伸")
            if control_fail:
                detail.append(control_fail)
            print(f"[FAIL] {name} ({tname}): " + "；".join(detail))
        else:
            print(f"[ OK ] {name} ({tname}): 唯一位点 {perfect[0][1]:,}，"
                  f"无 3′ 端可退火错配{control_note} | Tm {tm(bind):.1f} °C")

    pair_rows = []
    pair_failed = 0
    byname = {p.get("name"): p for p in spec["primers"]}
    for pr in spec.get("pairs", []):
        pname = pr.get("name", "?")
        fps = byname.get(pr.get("f"))
        rps = byname.get(pr.get("r"))
        if not fps or not rps:
            print(f"[SKIP] 引物对 {pname}: 找不到 {pr.get('f')} 或 {pr.get('r')}")
            continue
        fseq = fps["seq"].upper()[-int(fps.get("binding", len(fps["seq"]))):]
        rseq = rps["seq"].upper()[-int(rps.get("binding", len(rps["seq"]))):]
        tname = pr.get("template")
        if tname not in tmpl_paths:
            print(f"[SKIP] 引物对 {pname}: 没有模板 '{tname}'")
            continue
        if tname not in cache:
            t0 = read_fasta(tmpl_paths[tname])
            cache[tname] = (t0, build_index(t0), tname not in linear)
        t0, idx0, circ0 = cache[tname]
        sizes, _, _ = pair_products(t0, idx0, fseq, rseq, args.max_mismatch, circ0)
        absent = list(pr.get("absent_on") or [])
        if pr.get("control"):
            absent.append(pr["control"])
        absent = [x for x in dict.fromkeys(absent) if x != tname]
        expect = pr.get("expect")
        tdesc = "、".join(f"{x:,} bp" for x in sizes) if sizes else "无产物"
        verdict = []
        if len(sizes) > 1:
            verdict.append("目标模板上有多个产物")
        if not sizes:
            verdict.append("目标模板上无产物")
        if expect is not None and sizes and sizes[0] != int(expect):
            verdict.append(f"大小与预期 {int(expect):,} bp 不符")
        bad = []
        for cname in absent:
            if cname not in tmpl_paths:
                bad.append(f"{cname}（没有给出序列文件）")
                verdict.append(f"对照 {cname} 未提供，无法确认")
                continue
            if cname not in cache:
                c0 = read_fasta(tmpl_paths[cname])
                cache[cname] = (c0, build_index(c0), cname not in linear)
            ct0, cidx0, ccirc0 = cache[cname]
            csizes, _, _ = pair_products(ct0, cidx0, fseq, rseq, args.max_mismatch, ccirc0)
            if csizes:
                bad.append(f"{cname} " + "、".join(f"{x:,} bp" for x in csizes))
                verdict.append(f"对照 {cname} 也有产物，不能区分")
        ok = len(sizes) == 1 and not bad and (expect is None or sizes[0] == int(expect))
        if bad:
            cdesc = "；".join(bad) + " 有产物"
        elif absent:
            shown = "、".join(absent[:3]) + ("…" if len(absent) > 3 else "")
            cdesc = f"{shown} {'均' if len(absent) > 1 else ''}无产物"
        else:
            cdesc = "未指定"
        pair_rows.append(dict(name=pname, target=tdesc, control=cdesc,
                              expect=expect if expect is not None else "",
                              verdict="；".join(verdict) or "可区分"))
        if ok:
            print(f"[ OK ] 引物对 {pname}：{tname} {tdesc}；对照 {cdesc}")
        else:
            pair_failed += 1
            print(f"[FAIL] 引物对 {pname}：{tname} {tdesc}；对照 {cdesc}"
                  f" -> {'；'.join(verdict)}")

    if pair_rows:
        print(f"\n{len(pair_rows)} 对引物，{pair_failed} 对不合格")
    print(f"{len(rows)} 条引物，{failed} 条不合格")
    if args.csv:
        import csv as _csv
        with open(args.csv, "w", newline="", encoding="utf-8-sig") as fh:
            w = _csv.writer(fh)
            w.writerow(["引物", "模板", "长度", "完全匹配", "唯一位置", "风险位点", "Tm(°C)", "GC(%)"])
            for r in rows:
                w.writerow([r["name"], r["template"], r["length"], r["exact"], r["pos"],
                            r["other"], r["tm"], r["gc"]])
            for r in pair_rows:
                w.writerow([f"引物对 {r['name']}", f"目标 {r['target']}",
                            f"对照 {r['control']}", f"预期 {r['expect'] or '-'}",
                            r["verdict"], "", "", ""])
        print(f"CSV -> {args.csv}")
    return 1 if (failed or pair_failed) else 0


if __name__ == "__main__":
    sys.exit(main())
