#!/usr/bin/env python3
"""Generate the lab primer sheet (xlsx) in the lab's own template layout.

Why a dedicated generator: the deliverable is a sheet the bench actually opens, with
formulas so the operator can change a primer and see 长度 / tm / 延伸时间 update.
Hand-typing the same numbers into Excel is how tables drift from the design.

Template layout (assets/primer_design_template.xlsx), one column set:

    A 产物      B 模板      C 名称      D 序列 5'-3'
    E 长度      F Tm        G tm        H 产物长度     I 延伸时间

One primer pair per block of two rows. A, B, G, H, I are merged over the pair.
E = LEN(D), G = AVERAGE(F pair) - 5, I = IF(H<=10000, H/1000*15, H/1000*30).

Input JSON:

    {
      "out": "引物表.xlsx",
      "groups": [
        {"product": "V16", "template": "pN1-16",
         "f": "VEC-F", "r": "VEC-R-A", "len": 5790},
        ...
      ],
      "primers": {
        "VEC-F": {"seq": "ACTCTGGGGTTCGAAATGACCGA", "tm": 61.5},
        "VEC-R-A": {"seq": "GCGTTCGTACTGTTCCACGATGG", "tm": 61.8}
      }
    }

`len` is the PCR product length on that template. `tm` is the Tm of the 3' binding
region only (scripts/tm_engine.py), never the whole primer with its 5' tail.

Exit code is non-zero when a primer referenced by a group is missing from `primers`
or carries no sequence, so a half-filled sheet cannot be shipped by accident.
"""

from __future__ import annotations

import argparse
import json
import sys

from openpyxl import Workbook
from openpyxl.styles import Alignment, Border, Font, Side
from openpyxl.utils import get_column_letter

HEADER = ["产物", "模板", "名称", "序列 5'-3'", "长度", "Tm", "tm", "产物长度", "延伸时间"]
WIDTHS = {"A": 16.0, "B": 13.3, "C": 15.1, "D": 51.3, "E": 9.1, "F": 9.1}
FONT_NAME = "等线"


def build(design: dict, path: str) -> int:
    groups = design.get("groups") or []
    primers = design.get("primers") or {}
    if not groups:
        print("no groups in design", file=sys.stderr)
        return 1

    missing, blank = [], []
    for g in groups:
        for key in ("f", "r"):
            name = g.get(key)
            if name not in primers:
                missing.append(name)
            elif not (primers[name].get("seq") or "").strip():
                blank.append(name)
    if missing or blank:
        if missing:
            print("primers missing from design:", ", ".join(sorted(set(map(str, missing)))), file=sys.stderr)
        if blank:
            print("primers without sequence:", ", ".join(sorted(set(blank))), file=sys.stderr)
        return 1

    wb = Workbook()
    ws = wb.active
    ws.title = "Sheet1"

    thin = Side(style="thin", color="000000")
    border = Border(left=thin, right=thin, top=thin, bottom=thin)

    for col, title in enumerate(HEADER, start=1):
        c = ws.cell(row=1, column=col, value=title)
        c.font = Font(name=FONT_NAME, size=12, bold=True)
        c.border = border
    for col, width in WIDTHS.items():
        ws.column_dimensions[col].width = width

    row = 2
    for g in groups:
        f_name, r_name = g["f"], g["r"]
        first, second = row, row + 1
        product = g.get("product", "")
        template = g.get("template", "")
        product_len = g.get("len")

        for offset, name in ((0, f_name), (1, r_name)):
            r = first + offset
            info = primers[name]
            ws.cell(row=r, column=1, value=product if offset == 0 else None)
            ws.cell(row=r, column=2, value=template if offset == 0 else None)
            ws.cell(row=r, column=3, value=name)
            ws.cell(row=r, column=4, value=info["seq"])
            ws.cell(row=r, column=5, value=f"=LEN(D{r})")
            ws.cell(row=r, column=6, value=info.get("tm"))
            ws.cell(row=r, column=7, value=f"=AVERAGE(F{first}:F{second})-5" if offset == 0 else None)
            ws.cell(row=r, column=8, value=product_len if offset == 0 else None)
            ws.cell(
                row=r,
                column=9,
                value=f"=IF(H{first}<=10000,H{first}/1000*15,H{first}/1000*30)" if offset == 0 else None,
            )

        for col in ("A", "B", "G", "H", "I"):
            ws.merge_cells(f"{col}{first}:{col}{second}")

        for r in (first, second):
            for col in range(1, 10):
                c = ws.cell(row=r, column=col)
                c.font = Font(name=FONT_NAME, size=11)
                letter = get_column_letter(col)
                if letter == "E":
                    c.number_format = '0"bp"'
                elif letter in ("F", "G"):
                    c.number_format = '0"℃"'
                elif letter == "I":
                    c.number_format = '0"s"'
                if letter in ("A", "B", "G", "H", "I"):
                    c.alignment = Alignment(horizontal="center", vertical="center", wrap_text=True)
        row += 2

    wb.save(path)
    print(f"wrote {path}: {len(groups)} pairs, {row - 2} data rows")
    return 0


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design", help="JSON with groups and primers")
    args = ap.parse_args(argv)
    with open(args.design, encoding="utf-8") as fh:
        design = json.load(fh)
    return build(design, design.get("out", "引物表.xlsx"))


if __name__ == "__main__":
    sys.exit(main())
