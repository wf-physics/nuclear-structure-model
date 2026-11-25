# -*- coding: utf-8 -*-
"""
F1 (S_n 闭式) 扩样评估: 10 测试核 → 85 核 (每条 Z 链近稳定线 1 核 + 原 10 测试核 + 滴线近邻)
留出集纪律: 闭式系数源自模型网格采样 (A=40..240, 零 AME 借用), 测试核天然不参与选参.
判据: |偏差|<=2 MeV 过 / <=5 半 / 其余不过 (与 q_formula_v4 一致)
预期 (论文口径): A>=40 区 MAE 1.30-1.58, 零"不过"; A<40 为闭式外推区精度降级
"""
import sys, os, json
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from nucdata import parse_ame

HERE = os.path.dirname(os.path.abspath(__file__))
ame = parse_ame(os.path.join(HERE, '..', 'data', 'mass_1.mas20'))
AV, AS, AC, ASYM = 13.850813715082808, 1.5924382285401872, 0.7264962440025764, 19.50

def Bf(Z, N):
    A = Z + N
    return AV*A - AS*A**(2/3) - AC*Z**2/A**(1/3) - ASYM*(N-Z)**2/A

def Sn_exp(Z, N):
    v1, v2 = ame.get((Z, N)), ame.get((Z, N-1))
    if not (v1 and v2 and v1.get('BEA') and v2.get('BEA')):
        return None
    return v1['BEA']*v1['A']/1000. - v2['BEA']*v2['A']/1000.

targets = []
for Z in range(8, 83):
    best = None
    for N in range(int(Z*1.0), int(Z*1.6)+2):
        if (Z, N) in ame and ame[(Z, N)].get('BEA') and (Z, N-1) in ame and ame[(Z, N-1)].get('BEA'):
            nz = abs(Z - N)
            if best is None or nz < best[0]:
                best = (nz, N)
    if best:
        targets.append((Z, best[1]))
extra = [(20,28),(28,30),(38,50),(50,62),(82,126),(6,8),(8,10),(17,20),(19,21),(25,30)]
targets = sorted(set(targets + extra))

rows = []
for Z, N in targets:
    se = Sn_exp(Z, N)
    if se is None:
        continue
    sm = Bf(Z, N) - Bf(Z, N-1)
    dev = sm - se
    v = '过' if abs(dev) <= 2.0 else ('半' if abs(dev) <= 5.0 else '不过')
    rows.append(dict(nucl=f"{ame[(Z,N)]['El']}-{Z+N}", Z=Z, N=N,
                     Sn_m=float(sm), Sn_e=float(se), dev=float(dev), v=v))

devs = np.array([abs(r['dev']) for r in rows])
ok = sum(1 for r in rows if r['v'] == '过'); hf = sum(1 for r in rows if r['v'] == '半')
print(f"F1 扩样: n={len(rows)}  过={ok}  半={hf}  不过={len(rows)-ok-hf}  "
      f"MAE={devs.mean():.2f}  最大偏差={devs.max():.2f}")
for lo, hi, tag in [(0,40,'A<40'), (40,140,'40<=A<140'), (140,250,'A>=140')]:
    t = [r for r in rows if lo <= r['Z']+r['N'] < hi]
    if t:
        d = np.array([abs(r['dev']) for r in t])
        o = sum(1 for r in t if r['v'] == '过'); h = sum(1 for r in t if r['v'] == '半')
        print(f"  {tag:12s}: n={len(t)} 过={o} 半={h} 不过={len(t)-o-h} MAE={d.mean():.2f}")
out = os.path.join(HERE, '..', 'results', 'f1_expanded_85.json')
json.dump(rows, open(out, 'w'), indent=1, ensure_ascii=False)
print(f"存档: {out}")
