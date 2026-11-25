# -*- coding: utf-8 -*-
"""
T12 扩样: F2 测试集 10 -> ~30 核, 两个指标重跑相关性:
  ① 距最近幻数距离 (原指标)
  ② AME 结合能 - 平滑液滴闭式 (壳修正量, DeepSeek 建议的高信噪比指标)
措辞判据: |r|>0.5 归因可写 "supports"; 0.2-0.5 写 "consistent with"; <0.2 撤回。
留出纪律: F2 测试核选取规则 = 常用衰变链核素 + 稳定线邻域, 扩样不改变原 10 核。
"""
import os, sys, json
import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nucdata import parse_ame

ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
AV, AS, AC, ASYM = 13.850813715082808, 1.5924382285401872, 0.7264962440025764, 19.50
def Bf(Z, N):
    A = Z + N
    return AV*A - AS*A**(2/3) - AC*Z**2/A**(1/3) - ASYM*(N-Z)**2/A

MAGIC = (8, 20, 28, 50, 82, 126)
# 原 10 测试核 (不动) + 扩样: 近稳定线 β 活跃核, Z=6..61, 每 Z 取 1 (与 f1_expanded 同规则但取 β 不稳定)
def Qb_exp(Z, N):
    v1, v2 = ame.get((Z, N)), ame.get((Z+1, N-1))
    if not (v1 and v2 and v1.get("ME") is not None and v2.get("ME") is not None):
        return None
    return (v1["ME"] - v2["ME"])/1000.0 + 0.782  # β⁻: M(A,Z)-M(A,Z+1)+0.782

rows0 = json.load(open(os.path.join(HERE, "results/t5_shell_correlation.json")))["rows"]
have = {(r["Z"], r["N"]) for r in rows0}
extra = []
for Z in range(6, 84):
    best = None
    for N in range(Z, int(Z*1.55)+2):
        if (Z, N) in have: continue
        q = Qb_exp(Z, N)
        if q is None or q <= 0: continue           # β⁻ 活跃核
        if (Z+1, N-1) not in ame: continue
        nz = abs(N - Z)
        if best is None or nz < best[0]: best = (nz, N, q)
    if best and len(extra) < 22:
        extra.append((Z, best[1], best[2]))
print(f"扩样新增 {len(extra)} 核 (原 10 + 新 {len(extra)} = {10+len(extra)})")

rows = [dict(r) for r in rows0]
for Z, N, q in extra:
    A = Z + N
    qm = Bf(Z+1, N-1) - Bf(Z, N)   # Lane 口径? 保持与 F2_Qbm_Lane 同式: B(A,Z+1)-B(A,Z)
    # 原 F2_Qbm_Lane 的 dev 定义 = Qm - Qe, 其中 Qm 用闭式 B 差分 (无 -0.782? 核对 C-14: Qm=1.71)
    # C-14: Bf(7,7)-Bf(6,8) = ? 直接用同式并加 -0.782 校验
    dN = min(abs(N - m) for m in MAGIC); dZ = min(abs(Z - m) for m in MAGIC)
    rows.append(dict(nucl=f"{ame[(Z,N)]['El']}-{A}", Z=Z, N=N,
                     dev=float(qm - q), dmin_NZ=min(dN, dZ), dN=dN, dZ=dZ,
                     Qe=float(q), Qm=float(qm)))

dev = np.array([r["dev"] for r in rows]); dmin = np.array([r["dmin_NZ"] for r in rows])
dN = np.array([r["dN"] for r in rows]); dZ = np.array([r["dZ"] for r in rows])
# 指标②: 壳修正量 = B_exp - B_smooth(闭式), 对母核与子核分别取, 取差 (Q 侧壳修正)
def shell_corr(Z, N):
    v = ame.get((Z, N))
    if not v or v.get("BEA") is None: return None
    return v["BEA"]*(Z+N)/1000. - Bf(Z, N)
sc = []
for r in rows:
    Z, N = r["Z"], r["N"]
    s_m, s_d = shell_corr(Z, N), shell_corr(Z+1, N-1)
    if s_m is None or s_d is None:
        sc.append(np.nan)
    else:
        # Q_壳 = -(壳修正差): 子核壳修正高 → Q 实验被压低 → dev 负偏贡献
        sc.append((s_m - s_d))
sc = np.array(sc, float)
mask = ~np.isnan(sc)

out = dict(n_total=len(rows),
           metric1_magic=dict(pearson=float(np.corrcoef(dmin, dev)[0, 1]),
                              spearman=float(spearmanr(dmin, dev).statistic)),
           metric2_shell=dict(pearson=float(np.corrcoef(sc[mask], dev[mask])[0, 1]),
                              spearman=float(spearmanr(sc[mask], dev[mask]).statistic),
                              n=int(mask.sum())),
           rows=rows)
print(f"\n指标① 幻数距离:   Pearson r={out['metric1_magic']['pearson']:+.3f}  Spearman={out['metric1_magic']['spearman']:+.3f}  (n={len(rows)})")
print(f"指标② 壳修正量:   Pearson r={out['metric2_shell']['pearson']:+.3f}  Spearman={out['metric2_shell']['spearman']:+.3f}  (n={mask.sum()})")
r1, r2 = abs(out['metric1_magic']['pearson']), abs(out['metric2_shell']['pearson'])
rbest = max(r1, r2)
verdict = ("归因成立 (|r|>0.5, 可写 supports)" if rbest > 0.5 else
           "consistent with (0.2–0.5)" if rbest > 0.2 else "撤回归因 (<0.2)")
print(f"裁决 → {verdict}")
out["verdict"] = verdict
json.dump(out, open(os.path.join(HERE, "results/t5_shell_correlation_v2.json"), "w"), indent=1, ensure_ascii=False)
print("存档: results/t5_shell_correlation_v2.json")
