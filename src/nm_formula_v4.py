# -*- coding: utf-8 -*-
"""
核物质标定 v4: 金刚石型四面体键网络 (四键饱和的无限极限) vs FCC 对照
物理动机: v5.0 引擎 k_max=4 四键饱和 → 四配位网络 → 金刚石构型 (四面体=α簇团, 与
已验证的 α 簇团结果同源); FCC 12 配位是单粒子图像, 与键网络模型身份不符。
方法: E 对 (k,C) 线性 → 每组 (α,p,δ) 解 2x2 方程组 [E(0.16)=-15.8, E'(0.16)=0], 扫 K
"""
import numpy as np, itertools, json

HB2M = 197.32698**2/(2*938.92)
RHOS = np.linspace(0.08, 0.40, 321)
RHO0 = 0.16
NPTS = 150

def diamond_sphere(N, RHO):
    """金刚石晶格球形截断: FCC 格点 + (1/4,1/4,1/4)a 基元, 8 原子/晶胞"""
    a = (8.0/RHO)**(1/3)
    m = int(np.ceil(a)) + 2
    pts = []
    fcc = [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]
    basis = [(0,0,0),(.25,.25,.25)]
    for i in range(-m,m+1):
        for j in range(-m,m+1):
            for k in range(-m,m+1):
                for f in fcc:
                    for b in basis:
                        p = ((i+f[0]+b[0])*a,(j+f[1]+b[1])*a,(k+f[2]+b[2])*a)
                        pts.append(p)
    P = np.array(pts)
    r = np.linalg.norm(P, axis=1)
    o = np.argsort(r)
    P = P[o][:N]
    P = P - P.mean(axis=0)
    return P

def fcc_sphere(N, RHO):
    a = (4.0/RHO)**(1/3)
    m = int(np.ceil(((3*N/(4*np.pi*RHO))**(1/3))/(a/2))) + 3
    idx = np.array([(i,j,k) for i,j,k in itertools.product(range(-m,m+1),repeat=3)
                    if (i+j+k)%2==0], float)
    P = idx*(a/2.0)
    r = np.linalg.norm(P, axis=1)
    return P[np.argsort(r)][:N]

def sig(x):
    return 1.0/(1.0+np.exp(np.clip(x, -60, 60)))

def fermi_T(rho):
    return 0.6*HB2M*(3*np.pi**2*rho/2)**(2.0/3)

def lattice_curves(mk_lattice, rhos):
    """返回 La(rho) (delta=0.02 基准) 与 D0 缩放接口"""
    P0 = mk_lattice(NPTS, RHO0)
    D0 = np.linalg.norm(P0[:,None,:]-P0[None,:,:], axis=2)
    np.fill_diagonal(D0, 1e9)
    scale = (RHO0/rhos)**(1.0/3)
    La = np.zeros(len(rhos))
    for i, rho in enumerate(rhos):
        D = D0*scale[i]
        La[i] = 0.5*((2.2-D)*sig((D-2.2)/0.02)*sig((0.8-D)/0.02)).sum()/NPTS
    def Lr_of(p):
        out = np.zeros(len(rhos))
        for i, rho in enumerate(rhos):
            D = D0*scale[i]
            out[i] = 0.5*((0.5/np.maximum(D,1e-6))**p).sum()/NPTS
        return out
    def Lad_of(delta):
        out = np.zeros(len(rhos))
        for i, rho in enumerate(rhos):
            D = D0*scale[i]
            out[i] = 0.5*((2.2-D)*sig((D-2.2)/delta)*sig((0.8-D)/delta)).sum()/NPTS
        return out
    return La, Lr_of, Lad_of

T = np.array([fermi_T(r) for r in RHOS])

def solve(lad, g, i0):
    dLad = np.gradient(lad, RHOS); dg = np.gradient(g, RHOS); dT = np.gradient(T, RHOS)
    A = np.array([[-lad[i0], g[i0]], [-dLad[i0], dg[i0]]])
    b = np.array([-15.8 - T[i0], -dT[i0]])
    try:
        k, C = np.linalg.solve(A, b)
    except np.linalg.LinAlgError:
        return None
    if k <= 0 or C <= 0:
        return None
    d2Lad = np.gradient(dLad, RHOS); d2g = np.gradient(dg, RHOS); d2T = np.gradient(dT, RHOS)
    E2 = -k*d2Lad[i0] + C*d2g[i0] + d2T[i0]
    K = 9*RHO0**2*E2
    eas = -k*lad + C*g + T
    return k, C, K, RHOS[np.argmin(eas)], eas[i0]

if __name__ == "__main__":
    T_arr = T
    summary = {}
    for name, mk in [("diamond(四键网络)", diamond_sphere), ("FCC(对照)", fcc_sphere)]:
        La, Lr_of, Lad_of = lattice_curves(mk, RHOS)
        i0 = np.argmin(np.abs(RHOS-RHO0))
        results = []
        for alpha in [0.5, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0]:
            for p in [2, 3, 4, 5, 6, 8]:
                for delta in [0.02, 0.05, 0.10, 0.20, 0.40]:
                    lad = Lad_of(delta)
                    g = (RHOS/RHO0)**alpha * Lr_of(p)
                    sol = solve(lad, g, i0)
                    if sol is None:
                        continue
                    k, C, K, rmin, e0 = sol
                    # 额外要求: 全局极小必须在 0.16±0.02 (否则解的是局域方程但全局极小跑偏)
                    if abs(rmin - 0.16) > 0.02:
                        continue
                    results.append((abs(K-240), alpha, p, delta, k, C, K, rmin, e0))
        results.sort()
        summary[name] = results
        print("="*78)
        print(f"{name}: 通过全局极小校验的解 {len(results)} 组")
        print(f"{'ΔK':>7s} {'α':>5s} {'p':>3s} {'δ':>5s} {'k':>8s} {'C':>9s} {'K':>8s} {'ρmin':>7s} {'E/A':>7s}")
        for dK, alpha, p, delta, k, C, K, rmin, e0 in results[:10]:
            print(f"{dK:7.0f} {alpha:5.2f} {p:3d} {delta:5.2f} {k:8.2f} {C:9.2f} {K:8.1f} {rmin:7.3f} {e0:7.2f}")
    if summary.get("diamond(四键网络)"):
        dK, alpha, p, delta, k, C, K, rmin, e0 = summary["diamond(四键网络)"][0]
        json.dump(dict(lattice="diamond", alpha=alpha, p=p, delta=delta, k=float(k),
                       C=float(C), K=float(K), rho0=float(rmin), eA=float(e0)),
                  open('results/nm_formula_v4.json','w'), indent=1)
        print("\n存档: results/nm_formula_v4.json")
