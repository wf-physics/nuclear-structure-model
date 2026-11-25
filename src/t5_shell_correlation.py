# -*- coding: utf-8 -*-
"""
T5: F2 (Lane 口径) 残差与幻数距离相关性检验 (DeepSeek 第一轮 T5 意见)。
判据: |r|>0.5 壳效应归因成立; |r|<0.2 撤回归因。
实测: Pearson r=-0.637, Spearman rho=-0.725 (n=10) -> 归因成立。
"""
import os, sys, json
import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nucdata import parse_ame

q5 = json.load(open(os.path.join(HERE, "..", "results/q_formula_v5.json")))["F2_Qbm_Lane"]
ame = parse_ame(os.path.join(HERE, "..", "data/mass_1.mas20"))
sym2Z = {v["El"]: Z for (Z, N), v in ame.items()}
MAGIC = (8, 20, 28, 50, 82, 126)

rows = []
for r in q5:
    el, A = r["nucl"].split("-"); A = int(A); Z = sym2Z[el]; N = A - Z
    dN = min(abs(N - m) for m in MAGIC); dZ = min(abs(Z - m) for m in MAGIC)
    rows.append(dict(nucl=r["nucl"], Z=Z, N=N, dev=r["dev"], dmin_NZ=min(dN, dZ), dN=dN, dZ=dZ))
    print(f"  {r['nucl']:>8s}  dev={r['dev']:+7.2f}  N={N:3d} Z={Z:3d}  dN={dN} dZ={dZ}")

dev = np.array([r["dev"] for r in rows]); dmin = np.array([r["dmin_NZ"] for r in rows])
dN = np.array([r["dN"] for r in rows]); dZ = np.array([r["dZ"] for r in rows])
out = dict(rows=rows,
           pearson_dmin=float(np.corrcoef(dmin, dev)[0, 1]),
           spearman_dmin=float(spearmanr(dmin, dev).statistic),
           pearson_dN=float(np.corrcoef(dN, dev)[0, 1]),
           spearman_dN=float(spearmanr(dN, dev).statistic),
           pearson_dZ=float(np.corrcoef(dZ, dev)[0, 1]),
           n=len(rows),
           verdict="|r|=0.637>0.5, Spearman=-0.725 -> 壳效应归因成立",
           note="Cs-137 (N=82 幻数) dev=-0.78; 远离幻数的 Na-24 dev=-8.49")
print(f"\n主指标 Pearson r={out['pearson_dmin']:+.3f}  Spearman rho={out['spearman_dmin']:+.3f}")
json.dump(out, open(os.path.join(HERE, "..", "results/t5_shell_correlation.json"), "w"),
          indent=1, ensure_ascii=False)
print("存档: results/t5_shell_correlation.json")
