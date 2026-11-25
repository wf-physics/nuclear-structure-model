# -*- coding: utf-8 -*-
"""
大验证: H→Ca 全量核素 × 四项观测量 (逐核素对错清单 + 正确率汇总)
====================================================================
数据源(全部官方):
  结合能 B     : AME2020 (BEA)
  电荷半径 r_ch: IAEA charge_radii.csv (Angeli/Marinova 2013)
  电四极矩 Q   : IAEA elec_mom_recomm.csv (推荐值, barn→fm²)
  α 簇团       : N=Z 偶偶 4n 核 (⁴He→⁴⁰Ca 标准 α 链, 理想簇数 = Z)
判定纪律(按用户要求: 数值量按误差带, 不按相等判):
  B   : 相对误差 ≤15% 过 / 15-30% 半过 / >30% 不过
  r_ch: 模型/实验 比 ∈ [0.7, 1.3] 过
  Q_s : 符号一致 且 数量级 |log10(Qm/Qe)| ≤ 0.7 (5倍) 过; 符号对但幅度超 = 半过
  α   : 找到簇数 ≥ 50% 理想簇数 过
优化器: L-BFGS basin hopping (v4 修复版), 标定 C=12, r_eq=1.608 (r_min=1.4)
"""
import numpy as np, json, os, re, csv, sys, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
_src = open(os.path.join(HERE, "relax_v2_test.py")).read()
exec(_src.split('print("泡利排斥饱和势标定')[0])

C_REP, R_EQ = 12.0, 1.608

def relax_nuc(Z, N, hops=None, seed=13):
    from scipy.optimize import minimize
    A = Z + N
    if A == 1: return np.zeros((1,3)), 0.0
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    a0 = R_EQ/np.sqrt(2); L = int(np.ceil(2.2*A**(1/3)/a0))+2
    sites = []
    for i in range(-L, L+1):
        for j in range(-L, L+1):
            for k in range(-L, L+1):
                if (i+j+k) % 2: continue
                sites.append(np.array([i, j, k])*a0)
    sites.sort(key=lambda s: np.linalg.norm(s))
    x0 = np.array(sites[:A])
    if hops is None: hops = 30 if A <= 12 else (20 if A <= 24 else 12)

    def fun(xf):
        x = xf.reshape(A, 3)
        F, E, _, _ = forces_energy(x, q, 4, C_REP, R_EQ)
        return E, -F.flatten()

    def local_min(x):
        r = minimize(fun, x.flatten(), jac=True, method="L-BFGS-B",
                     options=dict(maxiter=350, ftol=1e-10, maxls=40))
        E, xa = r.fun, r.x
        r2 = minimize(fun, xa, jac=True, method="L-BFGS-B",
                      options=dict(maxiter=150, ftol=1e-12, maxls=40))
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
    ar = float(np.sqrt(max(lam[2], 1e-9)/max(lam[0], 1e-9)))
    pidx = np.where(q == 1)[0]; nidx = np.where(q == 0)[0]
    D = np.sqrt(((x[:, None]-x[None, :])**2).sum(-1)) + np.eye(A)*9
    cands = []
    for pp in itertools.combinations(pidx, 2):
        for nn in itertools.combinations(nidx, 2):
            ids = pp + nn
            edges = [D[i, j] for i, j in itertools.combinations(ids, 2)]
            if max(edges) <= 2.2 and min(edges) >= 0.7:
                cands.append(ids)
    cands.sort(key=lambda ids: max(D[i, j] for i, j in itertools.combinations(ids, 2))/
               min(D[i, j] for i, j in itertools.combinations(ids, 2)))
    used, ncl = set(), 0
    for ids in cands:
        if not (set(ids) & used): used |= set(ids); ncl += 1
    rch = float(np.sqrt(((xp-xp.mean(0))**2).sum(1).mean())) if len(xp) > 1 else 0.0
    rms = float(np.sqrt(((dm)**2).sum(1).mean()))
    return Q0, ar, ncl, rms, rch

# ================================================================ 实验数据
RCH = {}
for row in csv.DictReader(open(os.path.join(HERE, "data/charge_radii_iaea.csv"))):
    z = int(row["z"])
    if z == 0: continue
    try: RCH[(z, int(row["a"]))] = float(row["radius_val"])
    except: pass

QEXP = {}
for row in csv.DictReader(open(os.path.join(HERE, "data/elec_mom_iaea.csv"))):
    z = int(row["z"])
    try: A = int(row["n.n+n.z"])
    except: continue
    if row.get("energy [keV]", "").strip() not in ("0", ""): continue
    qs = row.get("electric quadrupole [b]", "").strip()
    m = re.match(r"([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)", qs)
    if not m: continue
    Qb = float(m.group(1))
    sp = row.get("spin", "").strip()
    ms = re.match(r"(\d+(?:/\d+)?)", sp)
    if not ms: continue
    s = ms.group(1)
    J = float(eval(s)) if "/" not in s else float(eval(s.replace("/", "/")))
    if J < 1: continue                     # J<1 光谱四极矩恒 0, 无法对照
    QEXP[(z, A)] = (Qb*100.0, J)           # barn → fm²

# 核素全集: Z≤20, AME 有结合能
NUCS = []
for (Z, N), v in sorted(ame.items(), key=lambda kv: (kv[0][0]+kv[0][1], kv[0][0])):
    if Z > 20 or Z < 1: continue
    if v.get("BEA") is None: continue
    B_exp = v["BEA"]*(Z+N)/1000.0
    if B_exp < 0.2: continue               # 非束缚态无意义
    NUCS.append((Z, N, B_exp))
print(f"核素全集: {len(NUCS)} 个 (Z≤20, AME 实测结合能)")
print(f"实验覆盖: r_ch {len([1 for Z,N,_ in NUCS if (Z,Z+N) in RCH])} 个, "
      f"Q {len([1 for Z,N,_ in NUCS if (Z,Z+N) in QEXP])} 个\n", flush=True)

def rotf(J): return J*(2*J-1)/((J+1)*(2*J+3))
def region(Z): return "轻(Z≤8)" if Z <= 8 else ("中(9-14)" if Z <= 14 else "重(15-20)")

# ================================================================ 主循环
rows = []
CKPT = os.path.join(HERE, "results/big_validation_ckpt.json")
done = {}
if os.path.exists(CKPT): done = json.load(open(CKPT))

for i, (Z, N, B_exp) in enumerate(NUCS):
    key = f"{Z},{N}"; A = Z+N
    if key in done:
        r = done[key]
    else:
        x, E = relax_nuc(Z, N)
        q = np.concatenate([np.ones(Z), np.zeros(N)])
        Q0, ar, ncl, rms, rch = geom_observables(x, q)
        r = dict(E=E, Q0=Q0, ar=ar, ncl=ncl, rms=rms, rch=rch)
        done[key] = r
        if i % 12 == 0:
            json.dump(done, open(CKPT, "w"))
    nm_ = f"{ame[(Z,N)]['El']}-{A}"
    B_m = -r["E"]
    rel = abs(B_m-B_exp)/B_exp
    vB = "过" if rel <= 0.15 else ("半" if rel <= 0.30 else "不过")
    re_ = RCH.get((Z, A))
    rch = r["rch"]
    vR = None
    if re_:
        ratio = rch/re_
        vR = "过" if 0.7 <= ratio <= 1.3 else "不过"
    qe = QEXP.get((Z, A))
    vQ = None; Q0m = r["Q0"]
    if qe:
        Q0e, J = qe
        Qsm = Q0m*rotf(J)
        sgn_ok = np.sign(Qsm) == np.sign(Q0e)
        mag_ok = abs(np.log10(max(abs(Qsm), 1e-9)/max(abs(Q0e), 1e-9))) <= 0.7
        vQ = "过" if (sgn_ok and mag_ok) else ("半" if sgn_ok else "不过")
    vA = None
    if Z == N and A % 4 == 0 and Z >= 2:
        vA = "过" if r["ncl"] >= 0.5*(A//4) else "不过"
    rows.append(dict(name=nm_, Z=Z, N=N, A=A, region=region(Z),
        B_model=round(B_m,2), B_exp=round(B_exp,2), B_rel=round(rel,4), B_v=vB,
        rch_m=round(rch,3), rch_e=re_, rch_ratio=round(rch/re_,3) if re_ else None, R_v=vR,
        Q0_model=round(Q0m,2), Q0_exp=round(qe[0],2) if qe else None,
        J=qe[1] if qe else None, Q_v=vQ,
        alpha_n=r["ncl"], alpha_ideal=(A//4) if vA is not None else None, A_v=vA,
        rms=round(r["rms"],2), axis=round(r["ar"],2)))
    if i % 12 == 0:
        print(f"[{i+1}/{len(NUCS)}] {nm_:8s} B={B_m:8.1f}/{B_exp:8.1f} ({rel*100:5.1f}%) {vB}"
              + (f"  r_ch={rch:.2f}/{re_:.2f} {vR}" if re_ else "")
              + (f"  Q{vQ}" if vQ else ""), flush=True)
json.dump(done, open(CKPT, "w"))
json.dump(rows, open(os.path.join(HERE, "results/big_validation_rows.json"), "w"),
          ensure_ascii=False, indent=1)

# ================================================================ 汇总
def acc(rows, key):
    t = [r for r in rows if r[key]]
    n_ok = sum(1 for r in t if r[key] == "过")
    n_half = sum(1 for r in t if r[key] == "半")
    return n_ok, n_half, len(t)

print("\n" + "="*70)
print("大验证正确率汇总 (H→Ca, %d 核素)" % len(rows))
print("="*70)
for key, lbl in [("B_v", "结合能 B"), ("R_v", "电荷半径 r_ch"), ("Q_v", "四极矩 Q_s"), ("A_v", "α 簇团")]:
    n_ok, n_half, t = acc(rows, key)
    if t == 0: continue
    print(f"  {lbl:14s}: 完全过 {n_ok}/{t} = {n_ok/t*100:5.1f}%   (含半过 {n_ok+n_half}/{t} = {(n_ok+n_half)/t*100:5.1f}%)")

print("\n--- 分区正确率 (结合能; 检验'错的重核居多'假设) ---")
for reg in ["轻(Z≤8)", "中(9-14)", "重(15-20)"]:
    t = [r for r in rows if r["region"] == reg]
    n_ok = sum(1 for r in t if r["B_v"] == "过")
    n_half = sum(1 for r in t if r["B_v"] == "半")
    mean_rel = np.mean([r["B_rel"] for r in t])*100
    nt = len(t)
    print(f"  {reg:8s} n={nt:3d}  过={n_ok/nt*100:5.1f}%  过+半={(n_ok+n_half)/nt*100:5.1f}%  平均相对误差={mean_rel:5.1f}%")

print("\n--- Q_s 逐项对错清单 ---")
for r in rows:
    if r["Q_v"]:
        print(f"  {r['name']:8s} J={r['J']}  Q0模型={r['Q0_model']:8.2f}  Q_exp={r['Q0_exp']:8.2f} fm²  → {r['Q_v']}")

print("\n--- α 簇团逐项 (N=Z 4n 链) ---")
for r in rows:
    if r["A_v"]:
        print(f"  {r['name']:8s} α簇 {r['alpha_n']}/{r['alpha_ideal']}  → {r['A_v']}")

# CSV 逐核素清单
with open(os.path.join(HERE, "results/big_validation_table.csv"), "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print("\n逐核素清单: results/big_validation_table.csv")
