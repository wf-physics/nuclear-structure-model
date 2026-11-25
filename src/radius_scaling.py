#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""§2.4 电荷半径标度因子诊断：中位比、标度后带内比例与剩余散布。

输入: results/big_validation_rows.json (389 核大验证存档)
输出: results/radius_scaling_result.json
论文引用: 三层分工框架预印本英文版 v1.1 §2.4。
"""
import json
import numpy as np
import os

HERE = os.path.dirname(os.path.abspath(__file__))
rows = json.load(open(os.path.join(HERE, "..", "results", "big_validation_rows.json")))
rs = [(r["name"], r["rch_m"], r["rch_e"]) for r in rows
      if isinstance(r.get("rch_m"), (int, float)) and isinstance(r.get("rch_e"), (int, float))
      and r["rch_e"] > 0 and r["rch_m"] > 0]
arr = np.array([m / e for _, m, e in rs])
med = float(np.median(arr))
# 论文发表口径的中位比 0.794（97 核，含 2 个无模型半径的氢同位素）作为标度因子
F_PAPER = 0.794
scaled = arr / F_PAPER


def frac(a, lo, hi):
    return float(np.mean((a >= lo) & (a <= hi))) * 100


out = {
    "n_model_pairs": len(rs),
    "median_ratio_archived": med,
    "scale_factor_used": F_PAPER,
    "inband_0707_130_before": frac(arr, 0.7, 1.3),
    "inband_085_115_before": frac(arr, 0.85, 1.15),
    "inband_085_115_after": frac(scaled, 0.85, 1.15),
    "inband_070_130_after": frac(scaled, 0.7, 1.3),
    "scaled_median": float(np.median(scaled)),
    "scaled_mean_abs_dev": float(np.mean(np.abs(scaled - 1))),
    "scaled_std": float(np.std(scaled)),
}
json.dump(out, open(os.path.join(HERE, "..", "results", "radius_scaling_result.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
