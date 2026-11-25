# -*- coding: utf-8 -*-
"""
结构引擎 v2.2 诊断扫描（几何观测量验证的前置条件检验）
==========================================================
背景: 几何验证(Q0/α簇/轴比)要求物理合理的基态构型。
修复链: ①链状亚稳态 → FCC 密堆初始化  ②预注册参数无配对平衡点(数学bug)
        ③FCC 最近邻距离公式错误  ④格点数不足
本脚本: 修复后扫描 (C, r_eq) 全网格, 检验三重约束是否可同时满足:
   (a) ⁴He B/A ≈ 7.07   (b) ⁴⁰Ca rms ≈ 3.5 fm   (c) B/A 平台 ≈ 8.5
结论: 无解 —— 单尺度 pairwise 势 + 四键饱和缺密度依赖饱和项(泡利压强)。
配对平衡条件: dV/dr=0 → r_min = r_eq − F_WELL/(2C), 原标定 C=2,r_eq=1.45 时 r_min<0 (不存在)。
"""
import numpy as np, os, json

HERE = os.path.dirname(os.path.abspath(__file__))
_src = open(os.path.join(HERE, "relax_v2_test.py")).read()
exec(_src.split('print("泡利排斥饱和势标定')[0])

def relax_fcc(Z, N, C_rep, r_eq, restarts=5, iters=2500, seed=11):
    A = Z + N
    q = np.concatenate([np.ones(Z), np.zeros(N)])
    a0 = r_eq/np.sqrt(2)                       # FCC nn 距离 = a0·sqrt(2) = r_eq
    L = int(np.ceil(2.0*A**(1/3)/a0)) + 2
    sites = []
    for i in range(-L, L+1):
        for j in range(-L, L+1):
            for k in range(-L, L+1):
                if (i+j+k) % 2: continue
                sites.append(np.array([i, j, k])*a0)
    sites.sort(key=lambda s: np.linalg.norm(s))
    assert len(sites) >= A, f"sites {len(sites)} < {A}"
    sites = np.array(sites[:A])                # 最内层 A 个格点 = 紧致初始化
    rng = np.random.default_rng(seed)
    best = (np.inf, None)
    for s in range(restarts):
        x = sites + rng.uniform(-0.15, 0.15, (A, 3))
        qq = q[rng.permutation(A)]
        v = np.zeros_like(x)
        for t in range(iters):
            mom = 0.5 - 0.35*t/iters
            F, E, _, _ = forces_energy(x, qq, 4, C_rep, r_eq)
            v = mom*v + 0.01*F
            vn = np.linalg.norm(v, axis=1, keepdims=True)
            v = np.where(vn > 0.25, v*0.25/vn, v)
            x = x + v
            if t > iters*0.75:
                F2, _, _, _ = forces_energy(x, qq, 4, C_rep, r_eq)
                x = x + 0.004*F2
        _, E, _, _ = forces_energy(x, qq, 4, C_rep, r_eq)
        if E < best[0]: best = (E, x)
    return best[1], best[0]

if __name__ == "__main__":
    B_AME = {}
    for k in [(1,1),(1,2),(2,1),(2,2),(6,6),(8,8),(16,16),(20,20),(20,28)]:
        v = ame.get(k)
        if v and v.get("BEA"): B_AME[k] = v["BEA"]*v["A"]/1000.0
    out = []
    print(f"{'C':>4s} {'r_eq':>5s} {'r_min':>6s} | 标定Σ | B/A: 4He 12C 16O 36Ar 40Ca 48Ca | rms: 4He 12C 16O 40Ca")
    print(f"{'实验':>16s} |      |      7.07 7.68 7.98 8.52 8.55 8.55 |      1.68 2.3  2.7  3.5")
    for C_rep in [8.0, 10.0, 12.0]:
        for r_eq in [1.60, 1.70, 1.80]:
            rmin = r_eq - 5.0/(2*C_rep)
            if not (1.15 <= rmin <= 1.65): continue
            errs = []
            for k in [(1,1),(1,2),(2,1),(2,2)]:
                _, E = relax_fcc(k[0], k[1], C_rep, r_eq, restarts=3, iters=1200)
                errs.append(abs(-E - B_AME.get(k, 0)))
            row = dict(C=C_rep, r_eq=r_eq, r_min=rmin, calib_err=sum(errs), BA={}, rms={})
            line = f"{C_rep:4.0f} {r_eq:5.2f} {rmin:6.3f} | {sum(errs):5.1f}"
            for k in [(2,2),(6,6),(8,8),(16,16),(20,20),(20,28)]:
                x, E = relax_fcc(k[0], k[1], C_rep, r_eq)
                ba = -E/(k[0]+k[1]); row["BA"][f"{k[0]}-{k[1]}"] = ba
                line += f" {ba:5.2f}"
            line += " |"
            for k in [(2,2),(6,6),(8,8),(20,20)]:
                x, E = relax_fcc(k[0], k[1], C_rep, r_eq)
                dm = x - x.mean(0)
                rms = float(np.sqrt(((dm)**2).sum(1).mean()))
                row["rms"][f"{k[0]}-{k[1]}"] = rms
                line += f" {rms:5.2f}"
            out.append(row)
            print(line, flush=True)
    json.dump(out, open(os.path.join(HERE, "results/engine_v22_scan.json"), "w"), indent=1)
    print("\n已存 results/engine_v22_scan.json")
