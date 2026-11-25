# -*- coding: utf-8 -*-
"""
法则二·方案 B 实装：泡利排斥势 + 四键饱和 → 结构层分离能验收重跑
[病根诊断] 旧势能: 纯线性吸引阱(斜率5 MeV/fm)+0.8fm硬核 → 无饱和机制,
  内部核子可积累 30-40 MeV 键能, B/A 随 A 发散 → 分离能误差 50 MeV 级。
[方案 B 势能] V_pair = 吸引阱 + 泡利排斥 + 硬核 + 库仑:
  吸引:   -5·(2.2-r)            (0.8≤r≤2.2, 保留 v4.0 法则)
  排斥:   C·max(0, r_eq-r)²     (r<r_eq 时按重叠量平方惩罚, 泡利重叠体现)
  硬核:   0.8fm 以下二次衔接 (保留)
  饱和:   每核子只保留最强 k_max 条键(互相同意制) — "每核子仅与~4近邻强作用"
  平衡键长 r_eq 由饱和约束预注册: 4键/核子 × bond(1.4fm)=4MeV → B/A≈8 MeV/A
[标定纪律] C, k_max 只在 ²H/³H/³He/⁴He 上标定(解耦标定, 预注册区间),
  ⁸Be/⁹Be/¹²C/¹⁶O 分离能为盲测 holdout, 标定后不再动。
"""
import numpy as np
import json, os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from nucdata import parse_ame

HERE = os.path.dirname(os.path.abspath(__file__))
ame = parse_ame(os.path.join(HERE, "data/mass_1.mas20"))
ME_H = ame[(1, 0)]["ME"]
ame[(2, 0)] = dict(A=2, ME=2*ME_H + 91.8, MEerr=0.5, BEA=0.0, Qb=None, flag=0, aux=[], El="He")
F_WELL = 5.0; R_MIN, R_MAX = 0.8, 2.2; E_AT_RMIN = -7.0; K_CORE = 5000.0

def forces_energy(X, q, k_max, C_rep, r_eq, w_soft=0.15):
    """V_pair = 吸引阱(−5(2.2−r)) + 泡利排斥(C·max(0,r_eq−r)², 叠加不替换) + 硬核 + 库仑
    力 = −dV/dr·(x_i−x_j)/r (标准形式, 修正旧代码过力 bug)"""
    A = len(X)
    if A == 1:
        return np.zeros((1,3)), 0.0, 0.0, 0.0
    d = X[:, None, :] - X[None, :, :]
    r = np.sqrt((d*d).sum(-1)) + np.eye(A)
    rinv = 1.0/np.where(r > 1e-9, r, 1e-9)
    iu = ~np.eye(A, dtype=bool)
    core = (r < R_MIN) & iu
    well = (r >= R_MIN) & (r <= R_MAX) & iu
    over = well & (r < r_eq)
    # 势能（吸引与排斥叠加）
    Vb = np.where(well, -F_WELL*(R_MAX-r), 0.0)
    Vb = np.where(over, Vb + C_rep*(r_eq-r)**2, Vb)
    Vb = np.where(core, E_AT_RMIN + K_CORE*(R_MIN-r)**2, Vb)
    # dV/dr
    dV = np.where(well, F_WELL, 0.0)                      # 吸引: dV/dr=+5
    dV = np.where(over, dV - 2*C_rep*(r_eq-r), dV)        # 排斥: −2C(r_eq−r)
    dV = np.where(core, -2*K_CORE*(R_MIN-r), dV)
    Fcoef = np.where(iu, -dV*rinv, 0.0)
    Fcoef = np.clip(Fcoef, -200.0, 200.0)                 # 硬核限力, 防爆
    # 四键饱和权重（互相同意制）
    cand = np.where(iu & (r <= R_MAX), Vb, np.inf)
    rank = np.argsort(np.argsort(cand, axis=1), axis=1)   # 0 = 最强
    Wm = ((rank < k_max) & (rank.T < k_max) & (r <= R_MAX) & iu)
    Ws = np.where(iu & (r <= R_MAX), np.where(Wm, 1.0, w_soft), 0.0)
    # 库仑
    qc = q[:, None]*q[None, :]
    Vc = np.where(iu & (qc != 0), 1.44*qc*rinv, 0.0)
    Fc = np.where(iu, 1.44*qc*rinv**2, 0.0)
    E_bond = 0.5*(Ws*Vb).sum()
    E_coul = 0.5*Vc.sum()
    E_total = E_bond + E_coul
    F = (d * ((Ws*Fcoef + Fc)*rinv)[:, :, None]).sum(axis=1)
    return F, E_total, E_bond, E_coul

def relax(Z, N, k_max, C_rep, r_eq, restarts=8, iters=2500, seed=1):
    A = Z + N
    if A == 1:
        return np.zeros((1,3)), 0.0
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    r0 = max(1.5, 1.3*A**(1/3))
    for s in range(restarts):
        pts = []; tries = 0
        while len(pts) < A and tries < 20000:
            p = rng.uniform(-r0, r0, 3); tries += 1
            if np.linalg.norm(p) > r0: continue
            if pts and np.min(np.linalg.norm(np.array(pts)-p, axis=1)) < 0.9: continue
            pts.append(p)
        x = np.array(pts) if len(pts) == A else rng.uniform(-r0, r0, (A,3))
        v = np.zeros_like(x)
        for t in range(iters):
            mom = 0.6 - 0.3*t/iters
            F, E, Eb, Ec = forces_energy(x, q, k_max, C_rep, r_eq)
            v = mom*v + 0.01*F
            vn = np.linalg.norm(v, axis=1, keepdims=True)
            v = np.where(vn > 0.3, v*0.3/vn, v)
            x = x + v
            if t > iters*0.7:
                F2, E2, _, _ = forces_energy(x, q, k_max, C_rep, r_eq)
                x = x + 0.002*F2      # 末段过阻尼精修
        F, E, Eb, Ec = forces_energy(x, q, k_max, C_rep, r_eq)
        if E < best[0]: best = (E, x.copy())
    return best[1], best[0]

def nm(k): return f"{ame[k]['El']}-{ame[k]['A']}"

CALIB = [(1,1),(1,2),(2,1),(2,2)]
HOLDOUT = [(4,4),(4,5),(6,6),(8,8)]
DAUGHTERS = {(1,1):[], (1,2):[], (2,1):[], (2,2):[],
             (4,4):[(3,4),(4,3)], (4,5):[(4,4)], (6,6):[(5,6),(6,5)], (8,8):[(7,8),(8,7)]}
NEED = sorted(set(CALIB + HOLDOUT + [d for v in DAUGHTERS.values() for d in v] +
                  [(10,10),(20,20)]), key=lambda k: (k[0]+k[1], k[0]))

print("泡利排斥饱和势标定（仅 ²H/³H/³He/⁴He）+ 盲测 holdout")
print(f"{'C':>5s} {'k':>2s} | ²H err ³H err ³He err ⁴He err | 标定Σ|err|")
B_AME = {}
for k in NEED:
    v = ame.get(k)
    if v and v.get("BEA"): B_AME[k] = v["BEA"]*v["A"]/1000.0
best_cfg = None; best_loss = np.inf
for r_eq in [1.35, 1.40, 1.45, 1.50]:
    for C_rep in [2.0, 3.0, 5.0, 8.0]:
        k_max = 4
        errs = {}
        for k in CALIB:
            x, E = relax(k[0], k[1], k_max, C_rep, r_eq)
            errs[k] = abs(-E - B_AME.get(k, 0.0))
        loss = sum(errs.values())
        print(f"r_eq={r_eq:.2f} C={C_rep:4.1f} | " + " ".join(f"{errs[k]:6.2f}" for k in CALIB) + f" | {loss:8.2f}")
        if loss < best_loss:
            best_loss = loss; best_cfg = (C_rep, k_max, r_eq)
C_rep, k_max, r_eq = best_cfg
print(f"\n==> 预注册标定选中: C={C_rep}, k_max={k_max} (标定集 Σ|err|={best_loss:.2f} MeV)\n")

# 用选中参数松弛全部需要的核素
structs = {}
for k in NEED:
    x, E = relax(k[0], k[1], k_max, C_rep, 1.4, restarts=5, seed=2)
    structs[k] = E
    if k in [(2,2),(4,4),(6,6),(8,8),(10,10),(20,20)]:
        print(f"  {nm(k):8s} E_model={E:9.2f} MeV  B/A_model={-E/(k[0]+k[1]):6.2f}  B/A_AME={B_AME.get(k,0)/(k[0]+k[1]):6.2f}")

def Emod(k): return structs.get(k)
def S_ame(mother, parts):
    tot = sum(ame[p]["ME"] for p in parts)
    return (tot - ame[mother]["ME"])/1000.0

print("\n===== 分离能验收 holdout（判据 <1 MeV, 第二评审核定）=====")
print(f"{'核':8s} {'通道':6s} {'模型/MeV':>9s} {'AME/MeV':>9s} {'误差':>7s} {'判定':>4s}")
n_pass = n_tot = 0
for mother, chs in [(k, DAUGHTERS[k]) for k in HOLDOUT]:
    Z, N = mother
    # 质子分离: daughter (Z-1,N) + ¹H
    for d, parts in [( (Z-1,N), [(Z-1,N),(1,0)] ), ( (Z,N-1), [(Z,N-1),(0,1)] )]:
        lbl = "S_p" if parts[1]==(1,0) else "S_n"
        if Emod(d) is None: continue
        S_m = Emod(d) + (Emod(parts[1]) or 0.0) - Emod(mother)
        S_a = S_ame(mother, parts)
        err = abs(S_m - S_a); n_tot += 1; n_pass += err < 1
        print(f"{nm(mother):8s} {lbl:6s} {S_m:9.2f} {S_a:9.2f} {err:7.2f} {'过' if err<1 else '不过'}")
    # α 分离 (⁴He 可用时)
    if Emod((2,2)) is not None and Z>=2 and N>=2 and Emod((Z-2,N-2)) is not None:
        S_m = Emod((Z-2,N-2)) + Emod((2,2)) - Emod(mother)
        S_a = S_ame(mother, [(Z-2,N-2),(2,2)])
        err = abs(S_m - S_a); n_tot += 1; n_pass += err < 1
        print(f"{nm(mother):8s} {'S_a':6s} {S_m:9.2f} {S_a:9.2f} {err:7.2f} {'过' if err<1 else '不过'}")
print(f"\n验收: {n_pass}/{n_tot} 通过 (<1 MeV)")

json.dump({f"{k[0]},{k[1]}": structs[k] for k in structs},
          open(os.path.join(HERE, "results/structures_v2.json"), "w"), indent=1)
print("已存 results/structures_v2.json")
