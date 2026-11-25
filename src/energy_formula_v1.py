# -*- coding: utf-8 -*-
"""
能量公式体系 v1.0 — 核物质标定（v5.0 法则结构冻结版）
========================================================
法则结构冻结项（一字不改，来自 v5.0 法则二）：
  V_att(r) = -k*(2.2-r)          线性键阱, 窗口 0.8<=r<=2.2   [k 为公式层标定常数]
  V_rep(r) = C*(r_c/r)^p         泡利排斥核, r->0 发散        [公式层实现]
公式层新增项（审计整合 v2.0 采纳）：
  F(rho) = 1/(1+(rho/rhoc)^n)    密度饱和因子                  [s, rhoc, n]
  T      = (3/5)*hbar^2*k_F^2/2m, k_F=(3*pi^2*rho/2)^(1/3)   局域费米动能
  Lane   b*(tau_i·tau_j)/2       同位旋键项                    [b 由 a_sym 反解]
验收目标（预注册，跑完冻结）：
  rho0 = 0.16+-0.01 fm^-3, E/A = -15.8+-0.5 MeV, K = 240+-60 MeV
"""
import numpy as np, itertools, json, sys
from scipy.optimize import differential_evolution

HB2M = 197.32698**2 / (2*938.92)          # hbar^2/2m = 20.73 MeV·fm^2
RHO_TARGET, EA_TARGET, K_TARGET = 0.16, -15.8, 240.0
RHOS = np.linspace(0.10, 0.30, 21)        # E/A(rho) 曲线网格
NPTS = 100                                # FCC 球点数（标定用）

def fcc_sphere(N, RHO):
    """FCC 球形截断, 最靠近中心的 N 个格点 (确定性)"""
    a = (4.0/RHO)**(1/3)
    m = int(np.ceil(((3*N/(4*np.pi*RHO))**(1/3))/(a/2))) + 3
    idx = np.array([(i,j,k) for i,j,k in itertools.product(range(-m,m+1),repeat=3)
                    if (i+j+k)%2==0], float)
    P = idx*(a/2.0)
    r = np.linalg.norm(P, axis=1)
    return P[np.argsort(r)][:N]

# ---- 预编译: 每个 rho 的距离矩阵 (几何只依赖 rho, 与参数无关) ----
DS = {}
for rho in RHOS:
    P = fcc_sphere(NPTS, rho)
    D = np.linalg.norm(P[:,None,:]-P[None,:,:], axis=2)
    np.fill_diagonal(D, 1e9)
    DS[rho] = D

def potential(D, k, C, rc, p):
    """冻结结构: 线性键阱(窗口硬截止) + 幂律泡利核"""
    att = -k*(2.2-D)*((D>=0.8)&(D<=2.2))
    rep = C*(rc/np.maximum(D,1e-6))**p
    return att + rep

def fermi_T(rho):
    return 0.6*HB2M*(3*np.pi**2*rho/2)**(2.0/3)

def curve(par):
    k, C, rc, p, s, rhoc, n = par
    F = 1.0/(1.0+(RHOS/rhoc)**n)
    eas = []
    for i, rho in enumerate(RHOS):
        V = potential(DS[rho], k, C, rc, p)
        E_pair = 0.5*V.sum()*s*F[i]/NPTS          # 均匀物质: 全局 F(rho)
        eas.append(E_pair + fermi_T(rho))
    return np.array(eas)

def conditions(par, rhos=RHOS):
    eas = curve(par)
    i = int(np.argmin(eas))
    if i == 0 or i == len(rhos)-1:
        return 0.0, -99.0, 1e5
    lo, mi, hi = rhos[i-1], rhos[i], rhos[i+1]
    rho0 = mi - 0.5*(lo-hi)*(eas[i-1]-eas[i])/(eas[i-1]-2*eas[i]+eas[i])
    c2 = np.polyfit(rhos[max(0,i-2):i+3], eas[max(0,i-2):i+3], 2)
    return rho0, eas[i], 9*rho0**2*2*c2[0]

def loss(par):
    k, C, rc, p, s, rhoc, n = par
    if not (5<=k<=150 and C>0 and 0.05<=rc<=1.5 and 1<=p<=12
            and 0.1<=s<=4 and 0.02<=rhoc<=0.6 and 0.2<=n<=4):
        return 1e6
    rho0, e0, K = conditions(par)
    return ((rho0-RHO_TARGET)/0.01)**2 + ((e0-EA_TARGET)/0.5)**2 + ((K-K_TARGET)/60)**2

if __name__ == "__main__":
    bounds = [(10,120), (10,600), (0.05,1.5), (1,12), (0.1,4), (0.02,0.6), (0.2,4)]
    r = differential_evolution(loss, bounds, seed=7, maxiter=120, popsize=24,
                               tol=1e-4, polish=True, workers=-1)
    k, C, rc, p, s, rhoc, n = r.x
    rho0, e0, K = conditions(r.x)
    print("="*70)
    print("核物质标定结果 (v5.0 结构冻结, 7 个公式层常数)")
    print("="*70)
    print(f"  loss = {r.fun:.2f}")
    print(f"  线性井斜率 k   = {k:8.2f} MeV/fm   (v5.0 原值 5.0)")
    print(f"  泡利核  C,rc,p = {C:7.1f} MeV, {rc:.3f} fm, p={p:.2f}")
    print(f"  F(rho)  s,rhoc,n = {s:.3f}, {rhoc:.4f} fm^-3, n={n:.3f}")
    print(f"  ---- 三条件 ----")
    print(f"  rho0 = {rho0:.4f}  (目标 {RHO_TARGET}+-0.01)   {'✓' if abs(rho0-0.16)<=0.011 else '✗'}")
    print(f"  E/A  = {e0:.2f}   (目标 {EA_TARGET}+-0.5)    {'✓' if abs(e0+15.8)<=0.51 else '✗'}")
    print(f"  K    = {K:7.1f}  (目标 {K_TARGET}+-60)    {'✓' if abs(K-240)<=61 else '✗'}")
    np.save('results/nm_formula_v1.npy', r.x)
    json.dump(dict(loss=float(r.fun), k=float(k), C=float(C), rc=float(rc), p=float(p),
                   s=float(s), rhoc=float(rhoc), n=float(n),
                   rho0=float(rho0), eA=float(e0), K=float(K)),
              open('results/nm_formula_v1.json','w'), indent=1)
    print("  存档: results/nm_formula_v1.{npy,json}")
