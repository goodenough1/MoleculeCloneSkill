#!/usr/bin/env python3
"""对账：引物表、设计数据、验证目录、文档必须说的是同一件事。

交付物不是单个文件。改一条引物会牵连引物表、方案表格、插图、GenBank、验证目录，
只要有一处还停在旧版本，看着都像"做完了"。这个脚本把它们逐条比对，任何漂移都报出来。

    python3 scripts/reconcile.py --sheet 引物表.xlsx \
        --design A16.json --design A17.json --design TB7.json \
        --verify 验证 --doc 无缝克隆方案.tex

检查四项：

1. 表内同名引物的序列必须一致（防止一行改了一行没改）。
2. 设计 JSON 里的引物序列必须与表一致（表是唯一真源）。
3. 表里写的产物长度必须等于 `验证/<产物>.dna` 的实际长度。
4. `--doc` 里 `名称 & \\seqsplit{序列}` 形式的行必须与表一致。

退出码非 0 表示有不一致；交付前必须清零。
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import snapgene_io


def read_sheet(path):
    from openpyxl import load_workbook

    ws = load_workbook(path).active
    rows = []
    for r in range(2, ws.max_row + 1):
        name = ws.cell(r, 3).value
        if not name:
            continue
        product = ws.cell(r, 1).value or (rows[-1]["product"] if rows else "")
        rows.append(dict(product=product, name=name.strip(),
                         seq=(ws.cell(r, 4).value or "").strip().upper(),
                         size=ws.cell(r, 8).value))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sheet", required=True)
    ap.add_argument("--design", action="append", default=[], help="design JSON, repeatable")
    ap.add_argument("--verify", help="folder holding <product>.dna")
    ap.add_argument("--doc", help="tex file whose primer rows must match the sheet")
    args = ap.parse_args(argv)

    rows = read_sheet(args.sheet)
    by_name = {}
    problems = []
    for x in rows:
        if x["name"] in by_name and by_name[x["name"]] != x["seq"]:
            problems.append("表内 %s 有两种序列" % x["name"])
        by_name[x["name"]] = x["seq"]

    for path in args.design:
        d = json.load(open(path))
        for p in d.get("primers", []):
            want = by_name.get(p["name"])
            if want and want != p["seq"].upper():
                problems.append("%s 里 %s 与表不一致" % (os.path.basename(path), p["name"]))

    if args.verify and os.path.isdir(args.verify):
        found = {}
        for f in glob.glob(os.path.join(args.verify, "**", "*.dna"), recursive=True):
            found.setdefault(os.path.splitext(os.path.basename(f))[0], f)
        for x in rows:
            if not x["size"] or x["product"].endswith("-CHK"):
                continue
            f = found.get(x["product"])
            if not f:
                problems.append("验证目录里没有 %s.dna" % x["product"])
                continue
            seq, _ = snapgene_io.parse(f)
            if len(seq) != int(x["size"]):
                problems.append("%s 表里 %s bp，验证文件 %d bp" % (x["product"], x["size"], len(seq)))

    if args.doc and os.path.exists(args.doc):
        tex = open(args.doc, encoding="utf-8").read()
        for x in rows:
            m = re.search(re.escape(x["name"]) + r"\s*&\s*\\seqsplit\{([^}]*)\}", tex)
            if not m:
                problems.append("文档里没有 %s 这一行" % x["name"])
            elif m.group(1).upper() != x["seq"]:
                problems.append("文档里 %s 与表不一致" % x["name"])

    for p in problems:
        print("不一致:", p)
    print("对账完成：%d 条引物，%d 处不一致" % (len(rows), len(problems)))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
