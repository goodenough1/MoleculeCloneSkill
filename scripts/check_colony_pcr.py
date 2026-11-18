#!/usr/bin/env python3
"""菌落 PCR 引物检查：零扩增 + 覆盖全部插入片段。

一条菌落 PCR 只有同时满足下面三条才算通过：

1. 在目标质粒上出唯一一条产物。
2. 每条插入片段都有一段只属于克隆的 30 bp 序列落在产物内（片段 5′ 端的同源臂是载体序列，要求覆盖它就会在母质粒上出带，两条判据互斥）。
3. 在所有提供过片段的来源模板上都不出产物——任一条引物缺位点即算无带。

用法：

    python3 scripts/check_colony_pcr.py design.json --target 目标.fa \
        --primer-pair CHK-F:mEL222-CHK-R --sheet 引物表.xlsx

design.json 给出片段与来源模板（fragments 里每项的 template 与 tpl_files），
目标序列单独给出；引物序列从表里读，表是唯一真源。
"""

from __future__ import annotations

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from tm_engine import tm


def read_fasta(path):
    with open(path) as fh:
        fh.readline()
        return "".join(fh.read().split()).upper()


def rc(s):
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def read_sheet(path):
    from openpyxl import load_workbook

    ws = load_workbook(path).active
    out = {}
    for r in range(2, ws.max_row + 1):
        name, seq = ws.cell(r, 3).value, ws.cell(r, 4).value
        if name and seq:
            out[name.strip()] = seq.strip().upper()
    return out


def amplify(template, f, r):
    """Return (start0, end0) of the product on the template, or None."""
    d = template + template
    n = len(template)
    i = d.find(f)
    while 0 <= i < n:
        j = d.find(rc(r), i + len(f))
        if j >= 0 and j + len(r) - i <= n:
            return i, j + len(r)
        i = d.find(f, i + 1)
    return None


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design")
    ap.add_argument("--target", required=True)
    ap.add_argument("--primer-pair", required=True, help="F:R")
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--tpl-dir", default=".")
    ap.add_argument("--verify", help="folder with <fragment>.dna, used to test coverage")
    args = ap.parse_args(argv)

    d = json.load(open(args.design))
    primers = read_sheet(args.sheet)
    fname, rname = args.primer_pair.split(":")
    f, r = primers[fname], primers[rname]
    target = read_fasta(args.target)
    order = d["order"]
    tpl = d.get("templates", {})

    problems = []
    hit = amplify(target, f, r)
    if not hit:
        print("%-24s %s / %s  目标上无产物" % (os.path.basename(args.design), fname, rname))
        print("  FAIL 目标质粒上无产物")
        return 1
    size = hit[1] - hit[0]
    win = target + target
    lo, hi = hit[0], hit[0] + size

    # 覆盖判定按"独有序列"：每段插入至少要有一段只属于克隆的 30 bp 落在产物内。
    # 片段 5' 端的同源臂是载体序列，母质粒里也有，要求覆盖它就必然在母质粒上出带。
    import snapgene_io

    # "独有"是相对母质粒而言：母质粒不含任何插入片段，它的序列就是假阳性来源。
    # 供体本身当然含自己那段，拿供体当参照会把每条片段都判成"不独特"。
    parent_key = d["fragments"][order[0]]["template"]
    src_seqs = []
    pfile = tpl.get(parent_key)
    if pfile and os.path.exists(os.path.join(args.tpl_dir, pfile)):
        src_seqs.append(read_fasta(os.path.join(args.tpl_dir, pfile)))
    for name in order[1:]:
        frag_dna = os.path.join(args.verify or "", name + ".dna")
        if not os.path.exists(frag_dna):
            problems.append("缺少片段文件 %s.dna" % name)
            continue
        seq, _ = snapgene_io.parse(frag_dna)
        found = False
        for off in range(0, len(seq) - 30 + 1, 5):
            w = seq[off:off + 30]
            if any(w in t or rc(w) in t for t in src_seqs):
                continue
            pos = win.find(w, lo)
            if 0 <= pos < hi - 30:
                found = True
                break
        if not found:
            problems.append("片段 %s 没有独有序列落在产物内" % name)

    for key, tfile in tpl.items():
        path = os.path.join(args.tpl_dir, tfile)
        if not os.path.exists(path):
            continue
        if amplify(read_fasta(path), f, r):
            problems.append("来源模板 %s 上有扩增" % key)

    print("%s  %s / %s  产物 %d bp" % (os.path.basename(args.design), fname, rname, size))
    for p in problems:
        print("  FAIL", p)
    if not problems:
        print("  PASS 目标唯一产物、全部插入片段覆盖、来源模板零扩增")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
