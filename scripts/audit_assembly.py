#!/usr/bin/env python3
"""Pre-flight audit for multi-fragment scarless assembly designs (In-Fusion / Gibson / SLIC).

Checks, in order:
  1. every fragment can be sliced from the target and the designed order reassembles
     the target base-for-base (circular plasmids supported);
  2. terminal homologies do not collide: no two fragment ends share the same 25-mer
     unless they are the intended 3'->5' neighbour pair, and no fragment has identical
     5' and 3' ends (self-circularisation);
  3. homology arms stay in the 15-21 bp window, their Tm sits as close to 50 °C as
     the window allows, and every arm is unique in the reaction: no two arms may be
     identical, within two mismatches, or nested, and a duplicate copy elsewhere in
     the target may not fall near any fragment end;
  4. each primer's 3' binding region is unique on its template (both strands, circular);
  5. Tm of each primer's binding region (SantaLucia unified NN) and the annealing
     temperature per fragment.

Exit code is non-zero as soon as anything FAILs, so it can gate the rest of a pipeline.

Design file (JSON):

    {
      "target": "target.fa",            # or "target_seq": "ACGT..."
      "order": ["A1", "B1", "C1"],
      "fragments": {
        "A1": {"product": [3079, 319], "template": "p16"},
        "B1": {"product": [294, 1055]},
        "X1": {"seq": "ACGT...", "template": "synthetic"}
      },
      "templates": {"p16": "p16.fa"},
      "primers": [{"name": "P1-A1-F", "seq": "AGC...", "template": "p16",
                    "fragment": "A1", "binding": 35}]
    }

`product` is a 1-based inclusive interval on the target (or 0-based half-open with
`"coord": "0-based-halfopen"`); end < start wraps the origin. Homology arms are
checked against a 15-21 bp window and a 50 °C Tm target; override with --arm-min,
--arm-max, --arm-tm, --arm-tm-tol, --arm-tm-hard.
`seq` is accepted for fragments that are not a substring of the target (synthetic inserts).
`templates` and `primers` are optional. `binding` is the number of 3' bases that
anneal to the template (default 20); `fragment` groups primers so the script can
recommend an annealing temperature per PCR.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def rc(s: str) -> str:
    return s.translate(COMP)[::-1]


def read_fasta(path: str) -> str:
    parts = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            parts.append(line.upper())
    if not parts:
        raise SystemExit(f"empty FASTA: {path}")
    return "".join(parts)


# ---------------------------------------------------------------- sequence ops
def sub(seq: str, a: int, b: int, coord: str = "1-based-inclusive") -> str:
    """Slice a circular sequence.

    coord="1-based-inclusive" (default, SnapGene-style): [a, b], b may wrap.
    coord="0-based-halfopen" (Python-style): [a, b), b may wrap.
    """
    n = len(seq)
    d = seq + seq
    if coord == "0-based-halfopen":
        return d[a:b] if b > a else d[a:b + n]
    if b >= a:
        return d[a - 1:b]
    return d[a - 1:b + n]


def overlap(a: str, b: str, lo: int = 15, hi: int = 60):
    """longest suffix of a equal to prefix of b."""
    for k in range(min(hi, len(a), len(b)), lo - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return None


def same_circle(a: str, b: str) -> bool:
    return len(a) == len(b) and a in (b + b)


def circular_count(seq: str, probe: str) -> int:
    """Occurrences of probe on the plus strand of a circular sequence."""
    n = len(seq)
    d = seq + seq
    return sum(1 for i in range(n) if d[i:i + len(probe)] == probe)


# ---------------------------------------------------------------- Tm engine
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from tm_engine import tm as tm_nn, gc, annealing  # noqa: E402
from check_primers import build_index, check as strict_primer_check, scan  # noqa: E402


# ---------------------------------------------------------------- checks
def primer_sites(template: str, probe: str):
    """Strict uniqueness: one perfect site, no extendable <=2-mismatch site."""
    idx = build_index(template)
    return strict_primer_check(template, probe, idx)


def slice_fragments(target: str, fragments: dict, coord: str):
    out = {}
    for name, spec in fragments.items():
        if "seq" in spec:
            out[name] = spec["seq"].upper()
        elif "product" in spec:
            a, b = spec["product"]
            out[name] = sub(target, a, b, coord)
        else:
            raise SystemExit(f"fragment {name}: need 'product' or 'seq'")
    return out


def assemble(order, frags, lo, hi):
    merged = [(order[0], frags[order[0]], None)]
    cur = frags[order[0]]
    for nxt in order[1:]:
        k = overlap(cur, frags[nxt], lo, hi)
        if not k:
            raise ValueError(f"no overlap between {merged[-1][0]} and {nxt}: "
                             f"{cur[-25:]} .. {frags[nxt][:25]}")
        cur = cur + frags[nxt][k:]
        merged.append((nxt, frags[nxt], k))
    close = overlap(cur, frags[order[0]], lo, hi)
    if not close:
        raise ValueError(f"cannot close the circle: {cur[-25:]} .. {frags[order[0]][:25]}")
    return cur[:-close], merged, close


def end_audit(order, frags, end_len):
    ends = {}
    for name, seq in frags.items():
        ends[(name, "5")] = seq[:end_len]
        ends[(name, "3")] = seq[-end_len:]
    allowed = set()
    for i, name in enumerate(order):
        nxt = order[(i + 1) % len(order)]
        allowed.add(frozenset(((name, "3"), (nxt, "5"))))
    problems = []
    seen = set()
    for a, ea in ends.items():
        for b, eb in ends.items():
            if a == b or ea != eb:
                continue
            key = frozenset((a, b))
            if key in allowed:
                continue
            if key in seen:
                continue
            seen.add(key)
            problems.append((a, b, ea))
    self_hits = [(n, frags[n][:end_len], frags[n][-end_len:])
                 for n in frags if frags[n][:end_len] == frags[n][-end_len:]]
    return problems, self_hits


def fragment_intervals(target, order, frags):
    """0-based start and length of each fragment on the circular target."""
    n = len(target)
    d = target + target
    out = {}
    for name in order:
        p = d.find(frags[name])
        out[name] = (p % n, len(frags[name])) if p >= 0 else None
    return out


def locate_fragment(pos1, intervals, n, end_len):
    """Which fragment contains a 1-based target position, and whether it is near an end."""
    p = (pos1 - 1) % n
    for name, iv in intervals.items():
        if iv is None:
            continue
        start, length = iv
        off = (p - start) % n
        if off < length:
            near = off < end_len or off >= length - end_len
            return name, ("end" if near else "internal")
    return None


# ---------------------------------------------------------------- reporting
class Report:
    def __init__(self):
        self.rows = []
        self.failed = False

    def section(self, title):
        self.rows.append(("", title))

    def ok(self, msg):
        self.rows.append(("PASS", msg))

    def fail(self, msg):
        self.failed = True
        self.rows.append(("FAIL", msg))

    def warn(self, msg):
        self.rows.append(("WARN", msg))

    def render(self):
        width = max(len(m) for tag, m in self.rows if tag) if self.rows else 40
        lines = []
        for tag, msg in self.rows:
            if not tag:
                lines.append("")
                lines.append(msg)
            else:
                lines.append(f"  [{tag}] {msg}")
        return "\n".join(lines)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("design")
    ap.add_argument("--end-len", type=int, default=25, help="terminal homology window compared (default 25)")
    ap.add_argument("--min-overlap", type=int, default=15)
    ap.add_argument("--max-overlap", type=int, default=60)
    ap.add_argument("--arm-min", type=int, default=15, help="shortest homology arm (default 15)")
    ap.add_argument("--arm-max", type=int, default=21, help="longest homology arm (default 21)")
    ap.add_argument("--arm-tm", type=float, default=50.0, help="target arm Tm in °C (default 50)")
    ap.add_argument("--arm-tm-tol", type=float, default=5.0, help="soft deviation from --arm-tm")
    ap.add_argument("--arm-tm-hard", type=float, default=10.0, help="hard deviation from --arm-tm")
    ap.add_argument("--binding", type=int, default=20, help="primer 3' bases used for Tm (default 20)")
    ap.add_argument("--max-mismatch", type=int, default=2, help="second-site mismatch tolerance (default 2)")
    ap.add_argument("--json", action="store_true")
    args = ap.parse_args(argv)

    spec = json.load(open(args.design))
    if "target_seq" in spec:
        target = spec["target_seq"].upper()
    else:
        target = read_fasta(spec["target"])

    coord = spec.get("coord", "1-based-inclusive")
    frags = slice_fragments(target, spec["fragments"], coord)
    order = spec["order"]
    missing = [n for n in order if n not in frags]
    if missing:
        raise SystemExit(f"order refers to unknown fragments: {missing}")
    unused = [n for n in frags if n not in order]
    rep = Report()
    rep.section(f"目标: {len(target):,} bp | 片段: {len(frags)} | 顺序: {' -> '.join(order)}")
    if unused:
        rep.warn(f"design contains fragments not used in 'order': {unused}")

    # 1. assembly
    rep.section("1. 组装复核")
    try:
        product, merged, close = assemble(order, frags, args.min_overlap, args.max_overlap)
    except ValueError as exc:
        rep.fail(str(exc))
        product = None
    if product is not None:
        chain = [f"{m[0]}({len(m[1]):,})" for m in merged]
        rep.ok("片段长度: " + " + ".join(chain))
        overlaps = [m[2] for m in merged[1:]] + [close]
        rep.ok("接头长度: " + ", ".join(str(k) for k in overlaps))
        if any(k < args.min_overlap for k in overlaps):
            rep.fail(f"有接头短于 {args.min_overlap} bp: {overlaps}")
        if len(product) != len(target) or not same_circle(product, target):
            rep.fail(f"组装产物与目标不一致（{len(product):,} vs {len(target):,} bp）")
        else:
            rep.ok(f"组装产物与目标逐碱基一致（{len(product):,} bp 环状）")

    # 2. terminal homology collisions
    rep.section("2. 末端同源碰撞")
    problems, self_hits = end_audit(order, frags, args.end_len)
    if problems:
        for a, b, seq in problems:
            rep.fail(f"末端重复: {a[0]}-{a[1]} 与 {b[0]}-{b[1]} 都是 {seq}")
        rep.fail("同源臂不唯一：该接头会有第二种配对方式，正确产物会被错误拼接/自环竞争")
    else:
        rep.ok(f"所有 {2 * len(frags)} 个末端两两比较，除设计好的相邻接头外无重复")
    for n, s5, s3 in self_hits:
        rep.fail(f"{n} 的 5′ 与 3′ 端相同（{s5}），会自环")
    if not self_hits:
        rep.ok("没有片段自身两端相同")

    # 3. homology arms: length, Tm, uniqueness in the reaction
    arms = []
    if product is not None:
        rep.section(f"3. 同源臂（{args.arm_min}–{args.arm_max} bp，Tm 尽量贴近 {args.arm_tm:.0f} °C）")
        junctions = []
        for i in range(1, len(merged)):
            junctions.append((f"{merged[i - 1][0]}→{merged[i][0]}", merged[i][0], merged[i][2]))
        junctions.append((f"{order[-1]}→{order[0]}", order[0], close))
        tgt_idx = build_index(target)
        intervals = fragment_intervals(target, order, frags)
        for label, fname, k in junctions:
            arm = frags[fname][:k]
            t_arm = tm_nn(arm)
            g_arm = gc(arm)
            notes = []
            hard = False
            warn = False
            if not (args.arm_min <= k <= args.arm_max):
                notes.append(f"长度 {k} bp 不在 {args.arm_min}–{args.arm_max}")
                hard = True
            dt = t_arm - args.arm_tm
            if abs(dt) > args.arm_tm_hard:
                notes.append(f"Tm 偏离 {dt:+.1f} °C")
                hard = True
            elif abs(dt) > args.arm_tm_tol:
                notes.append(f"Tm 偏离 {dt:+.1f} °C")
                warn = True
            hits = scan(target, arm, tgt_idx, 2)
            iv = intervals.get(fname)
            expected = iv[0] if iv else None
            def other(h):
                return expected is None or (h[1] - 1) % len(target) != expected
            perfect = [h for h in hits if h[0] == "+" and h[2] == 0 and other(h)]
            near = [h for h in hits if h[0] == "+" and h[2] > 0 and other(h)]
            minus = [h for h in hits if h[0] == "-"]
            for h in perfect:
                where = locate_fragment(h[1], intervals, len(target), args.end_len)
                if where and where[1] == "end":
                    notes.append(f"{where[0]} 末端另有相同拷贝（{h[1]:,} 位），会竞争配对")
                    hard = True
                elif where:
                    notes.append(f"{where[0]} 内部另有相同拷贝（{h[1]:,} 位，不影响末端配对）")
                    warn = True
                else:
                    notes.append(f"目标 {h[1]:,} 位另有相同拷贝")
                    warn = True
            for h in near:
                where = locate_fragment(h[1], intervals, len(target), args.end_len)
                tag = "会竞争配对" if where and where[1] == "end" else "近似位点"
                notes.append(f"目标 {h[1]:,} 位有 {h[2]} 错配{tag}")
                if where and where[1] == "end":
                    hard = True
                else:
                    warn = True
            if minus:
                notes.append(f"反向互补 {len(minus)} 处")
                warn = True
            # best length inside the window, measured on the same junction
            best = None
            upper = min(args.arm_max, len(frags[fname]))
            for kk in range(args.arm_min, upper + 1):
                tt = tm_nn(frags[fname][:kk])
                if best is None or abs(tt - args.arm_tm) < abs(best[1] - args.arm_tm):
                    best = (kk, tt)
            hint = ""
            if best:
                if best[0] != k:
                    delta = best[0] - k
                    hint = (f"；窗口内最接近的接头是 {best[0]} bp（Tm {best[1]:.1f}），"
                            f"把上一条片段的 3′ 端{'延' if delta > 0 else '缩'} {abs(delta)} bp 即可")
                if abs(best[1] - args.arm_tm) > args.arm_tm_tol:
                    hint += (f"；这个位置在 {args.arm_min}–{args.arm_max} bp 内都做不到 "
                             f"{args.arm_tm:.0f}±{args.arm_tm_tol:.0f} °C，建议把接头挪到旁边的序列")
            msg = (f"{label}: {k} bp, Tm {t_arm:.1f} °C, GC {g_arm:.0f}%"
                   + ("；" + "；".join(notes) if notes else "") + hint)
            if hard:
                rep.fail(msg)
            elif warn:
                rep.warn(msg)
            else:
                rep.ok(msg)
            arms.append((label, arm))

        dup = False
        for i in range(len(arms)):
            for j in range(i + 1, len(arms)):
                a, b = arms[i][1], arms[j][1]
                li, lj = arms[i][0], arms[j][0]
                if a == b:
                    rep.fail(f"{li} 与 {lj} 的同源臂完全相同（{a}），会出现第二种配对")
                    dup = True
                elif len(a) == len(b):
                    mm = sum(1 for x, y in zip(a, b) if x != y)
                    if mm <= 2:
                        rep.warn(f"{li} 与 {lj} 只差 {mm} 个碱基，可能交叉退火")
                elif len(a) < len(b) and a in b:
                    rep.warn(f"{li} 的接头完全包含在 {lj} 的接头里")
                elif len(b) < len(a) and b in a:
                    rep.warn(f"{lj} 的接头完全包含在 {li} 的接头里")
        if not dup:
            rep.ok(f"{len(arms)} 个接头两两不同")

    # 4/5. primers
    primers = spec.get("primers") or []
    templates = spec.get("templates") or {}
    if primers:
        rep.section("4. 引物 3′ 端唯一性与 Tm")
        tpl_cache = {}

        def tpl_seq(key):
            if key not in tpl_cache:
                tpl_cache[key] = read_fasta(templates[key]) if key in templates else None
            return tpl_cache[key]

        have_tpl = any(p.get("template") in templates for p in primers)
        strict_cache = {}
        groups = {}
        for p in primers:
            name = p.get("name", "?")
            seq = p["seq"].upper()
            t = tpl_seq(p.get("template")) if have_tpl else None
            nb = int(p["binding"]) if "binding" in p else None
            if nb is None:
                # 无 binding 时：能整条完美匹配就按全长查，否则按 3′ 20 nt 查
                if t is not None:
                    if id(t) not in strict_cache:
                        strict_cache[id(t)] = build_index(t)
                    _ok, _perf, _risk, _ = strict_primer_check(t, seq, strict_cache[id(t)])
                    nb = len(seq) if _perf else min(args.binding, len(seq))
                else:
                    nb = min(args.binding, len(seq))
            bind = seq[-nb:]
            if t is not None:
                if id(t) not in strict_cache:
                    strict_cache[id(t)] = build_index(t)
                ok, perfect, risky, _hits = strict_primer_check(t, bind, strict_cache[id(t)])
            tm_val = tm_nn(bind)
            tail = len(seq) - nb
            tstr = f"3′ {nb} nt, Tm {tm_val:.1f} °C, GC {gc(bind):.0f}%"
            if tail > 0:
                tstr += f", 5′ 尾 {tail} nt"
            if t is None:
                rep.ok(f"{name}: {tstr}")
            elif ok:
                rep.ok(f"{name}: 唯一位点 {perfect[0][1]:,}，无 ≤2 错配可延伸位点 | {tstr}")
            else:
                detail = []
                if len(perfect) != 1:
                    detail.append(f"完整匹配 {len(perfect)} 处")
                for h in risky[:3]:
                    detail.append(f"{h[0]} 链 {h[1]:,} 位 {h[2]} 错配 {h[3]}")
                rep.fail(f"{name}: " + "；".join(detail) + f"（模板 {p['template']}）| {tstr}")
            frag = p.get("fragment")
            if frag:
                groups.setdefault(frag, []).append((tm_val, tail))

        if groups:
            rep.section("5. 推荐退火温度（min Tm − 3，钳制 50–72 °C）")
            for frag, vals in groups.items():
                if len(vals) < 2:
                    continue
                tms = [v[0] for v in vals]
                ta = annealing(min(tms), max(tms))
                note = ""
                if any(t > 30 for _, t in vals):
                    note = "  ← 有 >30 nt 的 5′ 尾巴（AT 富集时再降 2–4 °C 或做梯度）"
                rep.ok(f"{frag}: Ta {ta} °C（Tm min {min(tms):.1f}）{note}")

    # summary
    rep.section("结论")
    if rep.failed:
        rep.fail("有检查未通过，先改设计再写文档/下单")
    else:
        rep.ok("全部通过")

    if args.json:
        print(json.dumps({"failed": rep.failed, "rows": rep.rows}, ensure_ascii=False, indent=2))
    else:
        print(rep.render())
    return 1 if rep.failed else 0


if __name__ == "__main__":
    sys.exit(main())
