# -*- coding: utf-8 -*-
"""
路径1 预研: 引擎级配对信号探测 (不改 v5.0 法则, 纯算法层)
关键洞察: 同 A 的偶偶核与奇奇核共享同一个金刚石格点骨架 (A 相同 → 同一组格点),
         库仑退火只有质子分布不同 → 格点壳噪声在差分中大部分相消.
         这是模型自有能量中配对信号唯一可能露头的通道.
检验: B_engine(偶偶) − B_engine(奇奇)  vs  实验同差 (≈ 液滴对称能差 + 2δ/√A)
判据: 若模型差与实验差的相关/差值 < 2 MeV 量级 → 引擎已长出配对结构 (DeepSeek 方向成立)
      若差值仍是噪声量级 (>5 MeV 散布) → 引擎无奇偶机制, 方案2 补丁是唯一路线
"""
import numpy as np, json, sys
sys.path.insert(0, '.')
from nucdata import parse_ame
from q_engine_pairing_core import binding  # 独立引擎核心 (复制自 v4, 见下)
AME = parse_ame('data/mass_1.mas20')

def Bexp(Z, N):
    v = AME.get((Z, N))
    return v['BEA']*v['A']/1000.0 if v and v.get('BEA') else None

CASES = []
for A in (48, 56, 88, 100, 120, 152, 200, 208):
    Z0 = A//2
    # 取 |N-Z| 最小的 ee 核 (Z 偶)
    for Z in range(Z0-3, Z0+4):
        N = A - Z
        if Z % 2 == 0 and N % 2 == 0 and Z >= 8 and N >= 8:
            Z0 = Z
            break
    ee = (Z0, A-Z0)
    oo = (Z0+1, A-Z0-1)   # 奇奇 (Z 奇 N 奇)
    CASES.append((A, ee, oo))

print("="*84)
print("引擎级同 A 差分: B(ee) − B(oo) vs 实验 (噪声相消通道)")
print("="*84)
rows = []
for A, ee, oo in CASES:
    Be_m = binding(*ee)
    Bo_m = binding(*oo)
    Be_e, Bo_e = Bexp(*ee), Bexp(*oo)
    if None in (Be_m, Bo_m, Be_e, Bo_e):
        print(f"  A={A}: 数据缺失, 跳过"); continue
    dm = Be_m - Bo_m
    de = Be_e - Bo_e
    rows.append(dict(A=A, ee=ee, oo=oo, dm=float(dm), de=float(de), dev=float(dm-de)))
    print(f"  A={A:3d} ee{ee} oo{oo}: 模型差={dm:+8.2f} 实验差={de:+8.2f} 偏差={dm-de:+8.2f}")

devs = np.array([r['dev'] for r in rows])
print()
print(f"偏差: MAE={abs(devs).mean():.2f}  RMS={np.sqrt((devs**2).mean()):.2f}  σ={devs.std():.2f}")
verdict = ('引擎差分与实验一致 (引擎已含奇偶机制雏形) → 路径1 可行'
           if abs(devs).mean() < 3.0 else
           '引擎差分仍是噪声量级 → 引擎无独立奇偶机制, 方案2 半经验补丁为 F2 唯一可行路线')
print(f"裁决 → **{verdict}**")
json.dump(rows, open('results/engine_pairing_probe.json', 'w'), indent=1, ensure_ascii=False)
print("存档: results/engine_pairing_probe.json")
