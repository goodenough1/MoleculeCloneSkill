# MolecularCloneSkill

多片段无缝克隆（In-Fusion / Gibson / SLIC）的设计与上线前审计 skill：同源臂长度与 Tm、臂在反应内的唯一性、引物唯一结合位点、DpnI 处理、菌落 PCR 阴性对照，以及重复同源碰撞这类会悄悄失败的坑。

## 内容

- `SKILL.md`：工作流、硬规则、筛选引物设计、湿实验执行顺序、交付前自检
- `references/pitfalls.md`：22 个真实坑，含现象、根因、修法
- `scripts/audit_assembly.py`：组装、末端同源、同源臂（15–21 bp / Tm ≈ 50 °C / 反应内唯一）、引物唯一性、退火温度审计
- `scripts/check_primers.py`：引物唯一性硬检查、筛选引物对的产物预测与多模板对照（`absent_on`、`linear`）
- `scripts/tm_engine.py`：固定参数的 Tm 与退火温度引擎

## 用法

```bash
python3 scripts/audit_assembly.py design.json
python3 scripts/check_primers.py primers.json --csv primer_qc.csv
python3 scripts/tm_engine.py primers.json --csv primers.csv --latex
```

脚本的退出码非 0 表示有检查未通过，可以直接用来卡住后续步骤。

## 安装

```bash
git clone git@github.com:goodenough1/MolecularCloneSkill.git \
  ~/.claude/skills/infusion-assembly-design
```

依赖 Python 3.8+。`tm_engine.py` 有 Biopython 时用 Biopython 的 `Tm_NN`，没有时用内置的纯 Python 实现，两者数值一致。
