# -*- coding: utf-8 -*-
"""
几何观测量验证 v2 —— 自包含: 修复优化器后的重标定 + Q0/α簇/轴比/电荷半径验证
================================================================================
修复史:
  v1 链伪影 → v3 FCC紧致初始化 → v4 L-BFGS basin hopping (v3发现手写动量动力学
  不守恒: 紧致初态被硬核削顶过冲炸飞, E=+0.3 非束缚; 用户指正"纯算法错误", 属实)
标定纪律: 原 (C=2, r_eq=1.45) 标定被动力学 bug 污染(平衡点 r_min=0.2fm 压硬核,
  全体过度塌缩), 故用修复后优化器在同一标定集 (²H/³H/³He/⁴He) 重标定,
  并把配对平衡约束 r_min = r_eq − F/(2C) ∈ [1.15, 1.6] 写进网格(数学约束非拟合)。
"""
import numpy as np, json, os, itertools, sys

HERE = os.path.dirname(os.path.abspath(__file__))
_src = open(os.path.join(HERE, "relax_v2_test.py")).read()
exec(_src.split('print("泡利排斥饱和势标定')[0])

F_W = 5.0
def pair_rmin(C, r_eq): return r_eq - F_W/(2*C)

def relax_nuc(Z, N, C_rep, r_eq, hops=None, seed=13):
    """L-BFGS 局部极小 + basin hopping"""
    from scipy.optimize import minimize
    A = Z + N
    if A == 1: return np.zeros((1,3)), 0.0
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    a0 = r_eq/np.sqrt(2); L = int(np.ceil(2.2*A**(1/3)/a0))+2
    sites = []
    for i in range(-L, L+1):
        for j in range(-L, L+1):
            for k in range(-L, L+1):
                if (i+j+k) % 2: continue
                sites.append(np.array([i, j, k])*a0)
    sites.sort(key=lambda s: np.linalg.norm(s))
    x0 = np.array(sites[:A])
    if hops is None: hops = 40 if A <= 16 else (25 if A <= 28 else 16)

    def fun(xf):
        x = xf.reshape(A, 3)
        F, E, _, _ = forces_energy(x, q, 4, C_rep, r_eq)
        return E, -F.flatten()

    def local_min(x):
        r = minimize(fun, x.flatten(), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=400, ftol=1e-10, maxls=40))
        E, xa = r.fun, r.x
        r2 = minimize(fun, xa, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=200, ftol=1e-12, maxls=40))
        if r2.fun < E: E, xa = r2.fun, r2.x
        return E, xa.reshape(A, 3)

    rng = np.random.default_rng(seed)
    E0, x0m = local_min(x0)
    best_E, best_x = E0, x0m.copy()
    for h in range(hops):
        xh = best_x + rng.normal(0, 0.45, (A, 3))
        Eh, xh = local_min(xh)
        if Eh < best_E - 1e-6: best_E, best_x = Eh, xh
    return best_x, best_E

def geom_observables(x, q):
    A = len(x)
    xp = x[q == 1]
    if len(xp) >= 1:
        d = xp - xp.mean(0)
        T = sum(np.outer(3*r, r) - np.dot(r, r)*np.eye(3) for r in d)
        Q0 = float(np.linalg.eigvalsh(T)[2])
    else:
        Q0 = 0.0
    dm = x - x.mean(0)
    lam = np.sort(np.linalg.eigvalsh(dm.T @ dm))
    axis_ratio = float(np.sqrt(max(lam[2], 1e-9)/max(lam[0], 1e-9)))
    # α 簇团: 2p2n 六边全在 [0.7, 2.2] (键力程含硬核接触; 正则度只描述)
    pidx = np.where(q == 1)[0]; nidx = np.where(q == 0)[0]
    D = np.sqrt(((x[:, None]-x[None, :])**2).sum(-1)) + np.eye(A)*9
    cands = []
    for pp in itertools.combinations(pidx, 2):
        for nn in itertools.combinations(nidx, 2):
            ids = pp + nn
            edges = [D[i, j] for i, j in itertools.combinations(ids, 2)]
            if max(edges) <= 2.2 and min(edges) >= 0.7:
                cands.append((max(edges)/min(edges), ids))
    cands.sort()
    used, ncl = set(), []
    for reg, ids in cands:
        if not (set(ids) & used): used |= set(ids); ncl.append((reg, ids))
    rch = float(np.sqrt(((xp-xp.mean(0))**2).sum(1).mean())) if len(xp) > 1 else 0.0
    rms = float(np.sqrt(((dm)**2).sum(1).mean()))
    return Q0, axis_ratio, ncl, rms, rch

# ================================================================ 重标定
CALIB = [(1,1),(1,2),(2,1),(2,2)]
B_AME = {}
for k in CALIB:
    v = ame.get(k)
    if v and v.get("BEA"): B_AME[k] = v["BEA"]*v["A"]/1000.0

print("重标定 (修复后优化器, 平衡点约束 r_min∈[1.15,1.6])")
grid = []
for C in [8.0, 12.0, 16.0, 20.0]:
    for rmin_t in [1.2, 1.3, 1.4, 1.5]:
        r_eq = rmin_t + F_W/(2*C)
        if not (1.15 <= pair_rmin(C, r_eq) <= 1.6): continue
        grid.append((C, r_eq))
best = None; loss_b = np.inf
for C, r_eq in grid:
    errs = []
    for k in CALIB:
        _, E = relax_nuc(k[0], k[1], C, r_eq, hops=12, seed=5)
        errs.append(abs(-E - B_AME.get(k, 0.0)))
    if sum(errs) < loss_b: loss_b = sum(errs); best = (C, r_eq)
    print(f"  C={C:4.0f} r_eq={r_eq:.3f} r_min={pair_rmin(C,r_eq):.3f} | "
          + " ".join(f"{e:5.2f}" for e in errs) + f" | Σ={sum(errs):6.2f}", flush=True)
C_REP, R_EQ = best
print(f"==> 选中: C={C_REP}, r_eq={R_EQ:.3f} (r_min={pair_rmin(*best):.3f} fm)  Σ={loss_b:.2f}\n", flush=True)

# ================================================================ 全量几何验证
ALPHA_CHAIN = [(2,2),(4,4),(6,6),(8,8),(10,10),(12,12),(14,14),(16,16),(18,18),(20,20)]
ODD_Q = [(3,3),(3,4),(4,5),(5,5),(5,6),(6,7),(7,7),(8,9),(10,11),(11,12),(12,13),(13,14),(17,18),(17,20),(19,20)]
ALL = sorted(set(ALPHA_CHAIN + ODD_Q), key=lambda k: (k[0]+k[1], k[0]))

EXP = {  # Stone 2016 核矩表 (Q_spectroscopic, mb, J)
 (3,3):(-40.0, 1.5), (3,4):(53.0, 1.5), (5,5):(847.3, 3.0), (5,6):(40.7, 1.5),
 (7,7):(20.4, 1.0), (8,9):(-25.6, 2.5), (10,11):(91.3, 1.5), (11,12):(104.2, 2.5),
 (12,13):(199.4, 2.5), (13,14):(146.6, 2.5), (17,18):(-81.65, 1.5), (17,20):(-68.1, 1.5),
 (19,20):(58.9, 1.5)}
def rotf(J): return J*(2*J-1)/((J+1)*(2*J+3))
EXP_R = {(2,2):1.681,(3,3):2.390,(4,5):2.455,(6,6):2.470,(8,8):2.697,(10,10):3.005,
         (20,20):3.477,(6,7):2.435,(6,8):2.503,(20,28):3.474}

res = {}
print(f"{'核':7s} {'E':>9s} {'B/A模型':>7s} {'B/A实验':>7s} {'Q0模型':>8s} {'Q0实验':>8s} {'符':>3s} {'轴比':>5s} {'α簇':>6s} {'r_ch模型/实验':>13s}")
for k in ALL:
    x, E = relax_nuc(k[0], k[1], C_REP, R_EQ)
    q = np.concatenate([np.ones(k[0]), np.zeros(k[1])])
    Q0, ar, ncl, rms, rch = geom_observables(x, q)
    A_ = k[0]+k[1]; nm_ = ame[k]["El"]+"-"+str(A_)
    v = ame.get(k); be = v["BEA"]*A_/1000.0 if v and v.get("BEA") else None
    einfo = EXP.get(k)
    eQ0s = f"{einfo[0]/10/rotf(einfo[1]):8.2f}" if einfo else "     --"
    sgn = ("✓" if einfo and np.sign(Q0) == np.sign(einfo[0]) else "✗") if einfo else " --"
    re_ = EXP_R.get(k)
    print(f"{nm_:7s} {E:9.2f} {-E/A_:7.2f} {(be/A_ if be else 0):7.2f} {Q0:8.2f} {eQ0s:>8s} {sgn:>3s} {ar:5.2f} {len(ncl):3d}/{A_//4:<2d} {rch:5.2f}/"
          f"{re_ if re_ else '--'}", flush=True)
    res[f"{k[0]},{k[1]}"] = dict(name=nm_, E=E, Q0=Q0, axis_ratio=ar, rms=rms, r_ch=rch,
        alpha_clusters=[[round(r,3), list(map(int,i))] for r, i in ncl],
        B_AME=be, Q_exp=einfo[0] if einfo else None, J=einfo[1] if einfo else None)

print("\n===== [G1] Q0 汇总 =====")
pairs = [(k, res[f"{k[0]},{k[1]}"]) for k in EXP if f"{k[0]},{k[1]}" in res]
ok_s = sum(np.sign(r["Q0"]) == np.sign(EXP[k][0]) for k, r in pairs)
from scipy.stats import spearmanr
qm = np.array([r["Q0"] for _, r in pairs]); qe = np.array([EXP[k][0] for k, _ in pairs])
rho, p = spearmanr(np.abs(qm), np.abs(qe))
print(f"  符号一致率: {ok_s}/{len(pairs)}   |Q| 秩相关 ρ={rho:.3f} (p={p:.3f})")

print("\n===== [G2] 电荷半径 =====")
ok = tot = 0
for k, r in res.items():
    kk = tuple(map(int, k.split(",")))
    if kk in EXP_R and kk[0] >= 2:
        e = EXP_R[kk]; tot += 1
        ratio = r["r_ch"]/e
        ok += 0.7 < ratio < 1.3
        print(f"  {r['name']:7s} r_ch={r['r_ch']:5.2f} 实验={e:5.2f} 比={ratio:5.2f}")
print(f"  0.7-1.3 倍区间内: {ok}/{tot}")

print("\n===== [G3] α 簇团 =====")
for k in ALPHA_CHAIN:
    r = res[f"{k[0]},{k[1]}"]
    print(f"  {r['name']:7s} α簇 {len(r['alpha_clusters'])}/{(k[0]+k[1])//4}  轴比={r['axis_ratio']:.2f}  rms={r['rms']:.2f}")

res["_calib"] = dict(C=C_REP, r_eq=R_EQ, r_min=pair_rmin(*best), calib_err=loss_b)
json.dump(res, open(os.path.join(HERE, "results/geometry_validation.json"), "w"),
          indent=1, ensure_ascii=False)
print("\n已存 results/geometry_validation.json")
