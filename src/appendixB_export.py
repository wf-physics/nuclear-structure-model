# -*- coding: utf-8 -*-
"""
附录 B: 逐核子参数表导出 (论文大纲前置工作 #4)
- 正文表: ⁴He / ¹⁶O / ⁴⁸Ca 完整逐核子表 (编号/类型/坐标/k_i/邻居/r_ij/V_ij/泡利标记)
- 补充材料: 全量逐核子长表 CSV (315 核素) + 每核素汇总 CSV
势与键判定与 relax_v2_test.forces_energy 逐字一致 (引擎口径)。
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
F_WELL, R_MIN, R_MAX = 5.0, 0.8, 2.2
C_REP, R_EQ = 12.0, 1.608          # 几何引擎标定 (§5 口径)
K_CORE = 100.0                      # 与 relax_v2_test 一致 (若不同以引擎值为准)
try:
    import relax_v2_test as _rv
    K_CORE = getattr(_rv, "K_CORE", K_CORE)
except Exception:
    pass
W_SOFT = 0.15

S = json.load(open(os.path.join(HERE, "results/structures.json")))
ame_keys = {tuple(int(t) for t in k.split(",")): k for k in S}

def per_nucleon_table(Z, N, override=None):
    v = (override or S).get(f"{Z},{N}") or S[f"{Z},{N}"]
    X = np.array(v["coords"], float)
    A = len(X)
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    d = X[:, None, :] - X[None, :, :]
    r = np.sqrt((d*d).sum(-1)) + np.eye(A)
    iu = ~np.eye(A, dtype=bool)
    core = (r < R_MIN) & iu
    well = (r >= R_MIN) & (r <= R_MAX) & iu
    over = well & (r < R_EQ)
    Vb = np.where(well, -F_WELL*(R_MAX-r), 0.0)
    Vb = np.where(over, Vb + C_REP*(R_EQ-r)**2, Vb)
    Vb = np.where(core, -7.0 + 5000.0*(R_MIN-r)**2, Vb)   # E_AT_RMIN + K_CORE(R_MIN-r)^2, 与引擎逐字一致
    cand = np.where(iu & (r <= R_MAX), Vb, np.inf)
    rank = np.argsort(np.argsort(cand, axis=1), axis=1)
    Wm = ((rank < 4) & (rank.T < 4) & (r <= R_MAX) & iu)   # 互相同意制主键
    Ws = np.where(iu & (r <= R_MAX), np.where(Wm, 1.0, W_SOFT), 0.0)
    Vpair = Ws * Vb
    k = Wm.sum(axis=1)
    ksoft = ((Ws > 0) & ~Wm).sum(axis=1)
    Vc = np.where(iu & (q[:, None]*q[None, :] != 0), 1.44*(q[:, None]*q[None, :])/np.where(r > 1e-9, r, 1e-9), 0.0)
    rows = []
    for i in range(A):
        nbs = [j for j in range(A) if Wm[i, j] or Ws[i, j]]
        for j in nbs:
            rows.append(dict(
                nucleon=i+1, type="p" if q[i] == 1 else "n",
                x=round(X[i, 0], 3), y=round(X[i, 1], 3), z=round(X[i, 2], 3),
                k_main=int(k[i]), k_soft=int(ksoft[i]),
                neighbor=j+1, nb_type="p" if q[j] == 1 else "n",
                r_ij=round(float(r[i, j]), 3),
                V_ij=round(float(Vpair[i, j]), 3),
                bond_kind="main" if Wm[i, j] else "soft",
                pauli_push=bool(over[i, j]), hard_core=bool(core[i, j]),
                V_coul=round(float(Vc[i, j]/2), 3)))
    meta = dict(A=A, E_model=v["E_model"], E_bond=v["E_bond"], E_coul=v["E_coul"],
                Wb_min=v.get("Wb_min"), nn_min=v.get("nn_min"), rms=v.get("rms_radius"),
                k_dist={str(int(x)): int((k == x).sum()) for x in range(5)})
    return rows, meta

# 正文三核素: 统一用 §5 几何引擎现场弛豫 (C=12, r_eq=1.608, seed=13), 保证表内 E 与逐键能量自洽
# (全量长表仍用 structures.json 存档坐标 = 判定层引擎口径, 表注声明两引擎差异)
import sys as _sys
from geometry_validation import relax_nuc

def ensure_structure(Z, N):
    """正文表专用: 几何引擎现场弛豫 (seed=13 主口径), 结果存入 S_MAIN, 不覆盖全量表用的 S"""
    x, E = relax_nuc(Z, N, C_REP, R_EQ, seed=13)
    from geometry_validation import geom_observables
    qv = np.concatenate([np.ones(Z), np.zeros(N)])
    Q0, ar, ncl, rms, rch = geom_observables(x, qv)
    d = x
    rmat = np.sqrt(((d[:, None, :] - d[None, :, :])**2).sum(-1)) + np.eye(Z+N)
    ium = ~np.eye(Z+N, dtype=bool)
    S_MAIN[f"{Z},{N}"] = dict(Z=Z, N=N, A=Z+N, El=("He" if Z == 2 else "O" if Z == 8 else "Ca"),
                              E_model=float(E), E_bond=None, E_coul=None, coords=x.tolist(),
                              Wb_min=None, Wb_max=None, EC_max=None,
                              nn_min=float(rmat[ium].min()), rms_radius=float(rms))
    print(f"  (几何引擎弛豫 Z={Z},N={N}: E={E:.2f}, B/A={-E/(Z+N):.2f})")

S_MAIN = {}
for Z, N in [(2, 2), (8, 8), (20, 28)]:
    ensure_structure(Z, N)
    rows, meta = per_nucleon_table(Z, N, override=S_MAIN)
    nm = S_MAIN[f"{Z},{N}"]["El"] + "-" + str(Z+N)
    out = dict(name=nm, meta=meta, table=rows)
    json.dump(out, open(os.path.join(HERE, f"results/appendixB_{nm}.json"), "w"),
              indent=1, ensure_ascii=False)
    print(f"{nm}: {len(rows)} 键行, k 分布={meta['k_dist']}, E={meta['E_model']:.2f}")

# 全量逐核子长表 + 汇总
import csv
allrows, summary = [], []
for k, v in sorted(S.items(), key=lambda kv: (kv[1]["A"], kv[1]["Z"])):
    Z, N = v["Z"], v["N"]
    if Z + N == 1:
        summary.append(dict(nucl=v["El"]+"-1", Z=Z, N=N, A=1, E_model=v["E_model"],
                            E_bond=0, E_coul=0, n_bond_rows=0, k_dist="{}",
                            Wb_min=None, nn_min=None, rms=None))
        continue
    rows, meta = per_nucleon_table(Z, N)
    for r in rows:
        r = dict(r); r["nucl"] = v["El"] + "-" + str(Z+N)
        allrows.append(r)
    summary.append(dict(nucl=v["El"]+"-"+str(Z+N), Z=Z, N=N, A=Z+N,
                        E_model=round(v["E_model"], 3) if v["E_model"] is not None else None,
                        E_bond=round(v["E_bond"], 3) if v.get("E_bond") is not None else None,
                        E_coul=round(v["E_coul"], 3) if v.get("E_coul") is not None else None,
                        n_bond_rows=len(rows),
                        k_dist=json.dumps(meta["k_dist"]), Wb_min=v.get("Wb_min"),
                        nn_min=v.get("nn_min"), rms=round(v.get("rms_radius") or 0, 3)))

with open(os.path.join(HERE, "results/appendixB_all_nucleons.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(allrows[0].keys()))
    w.writeheader(); w.writerows(allrows)
with open(os.path.join(HERE, "results/appendixA_summary.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
    w.writeheader(); w.writerows(summary)
print(f"\n全量逐核子表: {len(allrows)} 行 → appendixB_all_nucleons.csv")
print(f"汇总表: {len(summary)} 核素 → appendixA_summary.csv")
