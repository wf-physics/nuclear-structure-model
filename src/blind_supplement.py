# -*- coding: utf-8 -*-
"""补充对照: 完全不拟合(随机参数)下模型的准确率 + 分模式准确率分解"""
import numpy as np
from blind_test import (predict, acc, EXP, nN, RANGES, THETA0, POOL, ame)

rng = np.random.default_rng(0)
accs = []
for _ in range(2000):
    th = tuple(10**rng.uniform(*r) for r in RANGES)
    accs.append(acc(predict(th), np.arange(nN)))
accs = np.array(accs)
print(f"随机参数(不拟合) 2000 组: mean={accs.mean()*100:.1f}%  min={accs.min()*100:.1f}%  max={accs.max()*100:.1f}%")

p = predict(THETA0)
print("分模式准确率 (原始参数, 无补丁):")
for m in ["stable", "bm", "bp", "n", "p", "alpha"]:
    idx = np.where(EXP == m)[0]
    if len(idx):
        print(f"  实验={m:7s} n={len(idx):3d}  判对 {acc(p[idx], idx)*100:.1f}%")
# 判错明细（无补丁版, 原参数）
bad = np.where(p != EXP)[0]
print(f"\n无补丁判错 {len(bad)} 个:")
for i in bad:
    d = POOL[i]; k = d["key"]
    print(f"  {ame[k]['El']}-{ame[k]['A']}: 实验={EXP[i]}  模型={p[i]}")
