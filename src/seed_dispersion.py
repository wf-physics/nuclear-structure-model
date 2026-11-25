# -*- coding: utf-8 -*-
"""
退火种子散布实验 (DeepSeek T7-6): L-BFGS+basin hopping 的结果依赖初值与种子,
对代表性核素子集多种子重跑, 报告 B/A、α簇团计数、Q0、r_ch 的散布 → 写入 limitations 第 6 条。
引擎标定与 geometry_validation.py 主流程一致 (seed=5 标定选中后, 全量验证用默认 seed=13)。
"""
import os, sys, json
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
HERE = os.path.dirname(os.path.abspath(__file__))
from geometry_validation import relax_nuc, geom_observables

NUCS = [(2, 2), (6, 6), (10, 10), (14, 14), (20, 20), (28, 20)]
SEEDS = [13, 5, 42, 101, 777]   # 13 = 主流程默认种子(存档口径), 其余为散布测试
C, R_EQ = 12.0, 1.608           # geometry_validation 主流程选中值 (与存档一致)

rows = []
for Z, N in NUCS:
    per = dict(nucl=f"Z={Z},N={N}", seeds=[])
    for sd in SEEDS:
        x, E = relax_nuc(Z, N, C, R_EQ, seed=sd)
        q = np.concatenate([np.ones(Z), np.zeros(N)])
        Q0, ar, ncl, rms, rch = geom_observables(x, q)
        per["seeds"].append(dict(seed=sd, E=float(E), Q0=float(Q0),
                                 axis_ratio=float(ar), alpha_n=int(len(ncl)), r_ch=float(rch)))
    for key in ("E", "Q0", "axis_ratio", "alpha_n", "r_ch"):
        vals = np.array([s[key] for s in per["seeds"]], float)
        per[f"{key}_mean"] = float(vals.mean())
        per[f"{key}_std"] = float(vals.std())
    per["E_spread_pct"] = float(per["E_std"] / abs(per["E_mean"]) * 100)
    rows.append(per)
    print(f"Z={Z:2d},N={N:2d}: E散布±{per['E_std']:.2f} MeV ({per['E_spread_pct']:.1f}%)  "
          f"Q0 ±{per['Q0_std']:.2f}  α簇计数={sorted(s['alpha_n'] for s in per['seeds'])}  "
          f"r_ch ±{per['r_ch_std']:.3f} fm", flush=True)

allb = [r["E_spread_pct"] for r in rows]
print(f"\n汇总 (n={len(rows)}核 × {len(SEEDS)}种子): E 散布 mean={np.mean(allb):.1f}% max={np.max(allb):.1f}%")
for key in ("alpha_n", "Q0", "r_ch"):
    mx = max(rows, key=lambda r: r[f"{key}_std"])
    print(f"  {key}: max std = {mx[f'{key}_std']:.3f} ({mx['nucl']})")
json.dump(rows, open(os.path.join(HERE, "results/seed_dispersion.json"), "w"), indent=1, ensure_ascii=False)
print("存档: results/seed_dispersion.json")
