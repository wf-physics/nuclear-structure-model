# -*- coding: utf-8 -*-
"""
法则四原型：逐核素 (Z,N) 的结构坐标搜索与能量松弛
势能（v4.0 §1.1/§4.4 的可微化实现，已声明的简化）:
  - 原始键能势（未折扣 S(k)，软饱和只用于事后能量报告）:
      0.8 ≤ r ≤ 2.2:  V = -5·(2.2-r)  (即 E_bond(r)=7(2.2-r)/1.4，取负号势阱)
      r < 0.8:        V = -7 + 5000·(0.8-r)²  (硬核排斥，连续衔接)
      r > 2.2:        V = 0
  - 库仑: +1.44·q_i·q_j/r
松弛: 多起点随机初始化 + 带动量梯度下降（力场 -∇V），取最低能构型。
输出: results/structures.json（坐标、E_model、W_beta_min(a=0.8/1.0)、E_C_max 等）
"""
import numpy as np
import json, os, sys, time

sys.path.insert(0, os.path.dirname(__file__))
from nucdata import parse_ame

K_CORE = 5000.0     # 硬核二次势系数 (MeV/fm²)
F_WELL = 5.0        # 键能斜率: dE_bond/dr = 5 MeV/fm
R_MIN, R_MAX = 0.8, 2.2
E_AT_RMIN = -F_WELL * (R_MAX - R_MIN)   # = -7 MeV

DATA = os.path.join(os.path.dirname(__file__), "data")
OUT = os.path.join(os.path.dirname(__file__), "results", "structures.json")

def forces_energy(X, q):
    """返回 (F (A,3), E_total MeV, E_bond, E_coul)"""
    A = len(X)
    d = X[:, None, :] - X[None, :, :]           # r_ij = x_i - x_j
    r = np.sqrt((d * d).sum(-1)) + np.eye(A)    # 避免除零
    iu = ~np.eye(A, dtype=bool)
    # 键能势
    V = np.zeros((A, A))
    Fmag = np.zeros((A, A))
    in_well = (r >= R_MIN) & (r <= R_MAX) & iu
    in_core = (r < R_MIN) & iu
    V = np.where(in_well, -F_WELL * (R_MAX - r), V)
    Fmag = np.where(in_well, -F_WELL, Fmag)     # 吸引: 力沿 -r_ij 方向
    V = np.where(in_core, E_AT_RMIN + K_CORE * (R_MIN - r) ** 2, V)
    Fmag = np.where(in_core, 2 * K_CORE * (R_MIN - r), Fmag)  # 排斥
    # 库仑
    qc = q[:, None] * q[None, :]
    Vc = np.where(iu & (qc != 0), 1.44 * qc / r, 0.0)
    Fc = np.where(iu, 1.44 * qc / r ** 3, 0.0)  # F·(x_i-x_j) 形式的系数
    FmagC = Fc
    E_total = 0.5 * ((V * iu).sum() + Vc.sum())
    E_bond = 0.5 * (V * iu).sum()
    E_coul = 0.5 * Vc.sum()
    F = (d * (Fmag + FmagC)[:, :, None]).sum(axis=1)  # Σ_j coef·(x_i-x_j)
    return F, E_total, E_bond, E_coul

def relax(Z, N, restarts=6, iters=4000, dt=0.005, seed=0):
    """多起点松弛（致密初始化 + 冷却动量 + 末段过阻尼精修），返回最优构型"""
    A = Z + N
    rng = np.random.default_rng(seed)
    if A == 1:
        return dict(coords=[[0.0, 0.0, 0.0]], E_total=0.0, E_bond=0.0, E_coul=0.0)
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    best = None
    r0 = max(1.4, 1.05 * A ** (1 / 3))
    for s in range(restarts):
        # 致密初始化 + 拒绝采样（对距 >= 0.85 fm）
        pts = []
        tries = 0
        while len(pts) < A and tries < 30000:
            p = rng.uniform(-r0, r0, 3)
            tries += 1
            if np.linalg.norm(p) > r0:
                continue
            if pts and np.min(np.linalg.norm(np.array(pts) - p, axis=1)) < 0.85:
                continue
            pts.append(p)
        x = np.array(pts) if len(pts) == A else rng.uniform(-r0, r0, (A, 3))
        v = np.zeros_like(x)
        for t in range(iters):
            mom = 0.9 - 0.4 * t / iters          # 0.9 → 0.5 冷却
            F, E, Eb, Ec = forces_energy(x, q)
            v = mom * v + dt * F
            vn = np.linalg.norm(v, axis=1, keepdims=True)
            v = np.where(vn > 0.04, v * (0.04 / np.maximum(vn, 1e-12)), v)
            x = x + v
            x = x - x.mean(0)
        # 末段过阻尼精修（纯梯度下降，消余振）
        for t in range(600):
            F, E, Eb, Ec = forces_energy(x, q)
            x = x + np.clip(dt * 0.5 * F, -0.02, 0.02)
            x = x - x.mean(0)
        F, E, Eb, Ec = forces_energy(x, q)
        if best is None or E < best["E_total"]:
            best = dict(coords=x.round(3).tolist(), E_total=round(E, 4),
                        E_bond=round(Eb, 4), E_coul=round(Ec, 4))
    return best

def observables(x, Z, N, a_list=(0.8, 1.0)):
    """从构型提取法则二可观测量"""
    x = np.asarray(x)
    xp, xn = x[:Z], x[Z:]
    out = {}
    # W_beta: 每个中子的 (1/a + Σ_p 1/r)·1.44，取最小
    if len(xn) > 0 and len(xp) > 0:
        d = np.sqrt(((xn[:, None, :] - xp[None, :, :]) ** 2).sum(-1))
        sum_inv = (1.0 / d).sum(1)
        out["Wb_self_only"] = 1.44 / min(a_list)
        out["Wb_min"] = {str(a): round(1.44 * (1.0 / a + sum_inv.min()), 4) for a in a_list}
        out["Wb_max"] = {str(a): round(1.44 * (1.0 / a + sum_inv.max()), 4) for a in a_list}
    elif len(xn) > 0:
        out["Wb_min"] = {str(a): round(1.44 / a, 4) for a in a_list}
        out["Wb_max"] = out["Wb_min"]
    else:
        out["Wb_min"] = None
    # E_C: 每个质子的 Σ_q≠p 1.44/r，取最大
    if len(xp) > 0:
        if len(xp) == 1:
            ec = np.array([0.0])
        else:
            dp = np.sqrt(((xp[:, None, :] - xp[None, :, :]) ** 2).sum(-1))
            iu = ~np.eye(len(xp), dtype=bool)
            ec = (np.where(iu, 1.44 / dp, 0.0)).sum(1)
        out["EC_max"] = round(float(ec.max()), 4)
    else:
        out["EC_max"] = 0.0
    # 构型几何: 最近邻距离, rms 半径
    if len(x) > 1:
        d = np.sqrt(((x[:, None, :] - x[None, :, :]) ** 2).sum(-1)) + np.eye(len(x))
        iu = ~np.eye(len(x), dtype=bool)
        out["nn_min"] = round(float(d[iu].min()), 3)
    else:
        out["nn_min"] = None
    com = x.mean(0)
    out["rms_radius"] = round(float(np.sqrt(((x - com) ** 2).sum(1).mean())), 3)
    return out

def main():
    ame = parse_ame(os.path.join(DATA, "mass_1.mas20"))
    targets = sorted([(z, n) for (z, n) in ame if ame[(z, n)]["A"] <= 40])
    print(f"待松弛核素: {len(targets)} (A<=40)")
    t0 = time.time()
    results = {}
    for idx, (z, n) in enumerate(targets):
        b = relax(z, n, seed=1000 * z + n)
        obs = observables(b["coords"], z, n)
        rec = dict(Z=z, N=n, A=z + n, El=ame[(z, n)]["El"],
                   E_model=b["E_total"], E_bond=b["E_bond"], E_coul=b["E_coul"],
                   coords=b["coords"], **obs)
        results[f"{z},{n}"] = rec
        if (idx + 1) % 20 == 0:
            el = time.time() - t0
            print(f"[{idx+1}/{len(targets)}] {ame[(z,n)]['El']}-{z+n} "
                  f"E={b['E_total']:.2f} MeV  已用 {el:.0f}s", flush=True)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, ensure_ascii=False)
    print(f"\n完成: {len(results)} 核素, 总用时 {time.time()-t0:.0f}s → {OUT}")
    # 抽查几个著名核素
    for key in ["1,1", "2,2", "6,6", "20,20"]:
        r = results.get(key)
        if r:
            print(f"{r['El']}-{r['A']}: E_model={r['E_model']:.2f} MeV, "
                  f"nn_min={r['nn_min']} fm, rms={r['rms_radius']} fm, "
                  f"Wb_min={r.get('Wb_min')}, EC_max={r.get('EC_max')}")

if __name__ == "__main__":
    main()
