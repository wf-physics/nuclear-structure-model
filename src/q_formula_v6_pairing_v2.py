# -*- coding: utf-8 -*-
"""
配对项探测 v2 — 修复: (1) 原网格全偶偶核, 扩充奇奇核+奇A核采样;
(2) 独立引擎 (不 import q_formula_v4, 避免触发其全流程); (3) 结果存盘供离线分析.
判定: 三奇偶组的残差均值分离度 vs 噪声 → 模型中是否存在可提取的配对信号
"""
import numpy as np, json, sys
sys.path.insert(0, '.')
from nucdata import parse_ame
AME = parse_ame('data/mass_1.mas20')

CAL = json.load(open('results/nm_formula_v4.json'))
k, C, alpha, p, delta = CAL['k'], CAL['C'], CAL['alpha'], CAL['p'], CAL['delta']
RHO0 = 0.16
HB2M = 197.32698**2/(2*938.92)
E2 = 1.44

def sig(x):
    return 1.0/(1.0+np.exp(np.clip(x, -60, 60)))

def diamond_sites_smooth(A, rho=RHO0):
    a = (8.0/rho)**(1/3)
    m = int(np.ceil(a)) + 2
    R_c = (3*A/(4*np.pi*rho))**(1/3)
    fcc = [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]
    basis = [(0,0,0),(.25,.25,.25)]
    def pts_within(fac):
        pts = []
        for i in range(-m, m+1):
            for j in range(-m, m+1):
                for kk in range(-m, m+1):
                    for f in fcc:
                        for b in basis:
                            pt = ((i+f[0]+b[0])*a, (j+f[1]+b[1])*a, (kk+f[2]+b[2])*a)
                            if np.linalg.norm(pt) <= R_c*fac:
                                pts.append(pt)
        return np.array(pts)
    P = pts_within(1.02)
    if len(P) < A:
        for fac in (1.08, 1.15, 1.25, 1.40):
            P = pts_within(fac)
            if len(P) >= A:
                break
        if len(P) < A:
            return None
    r = np.linalg.norm(P, axis=1)
    P = P[np.argsort(r)][:A]
    if len(P) != A:
        P = P*(A/len(P))**(1/3)
    return P - P.mean(axis=0)

def binding(Z, N):
    A = Z + N
    P = diamond_sites_smooth(A)
    if P is None:
        return None
    D = np.linalg.norm(P[:, None, :]-P[None, :, :], axis=2)
    np.fill_diagonal(D, 1e9)
    att = -k*(2.2-D)*sig((D-2.2)/delta)*sig((0.8-D)/delta)
    rep = C*(0.5/np.maximum(D, 1e-6))**p
    E_pair = 0.5*(att+rep).sum()
    rho_p, rho_n = Z/A*RHO0, N/A*RHO0
    kp, kn = (3*np.pi**2*rho_p)**(1/3), (3*np.pi**2*rho_n)**(1/3)
    T = 0.6*HB2M*(Z*kp**2 + N*kn**2)
    bestEc = np.inf
    for seed in (42, 43, 44):
        rng = np.random.default_rng(seed)
        isP = np.zeros(A, bool); isP[rng.choice(A, size=Z, replace=False)] = True
        def coulomb(mask):
            pp = np.outer(mask, mask)
            np.fill_diagonal(pp, False)
            return 0.5*(E2/D[pp]).sum() if pp.any() else 0.0
        Ec = coulomb(isP)
        for it in range(6000):
            temp = 0.08*np.exp(-3.5*it/6000) + 0.01
            ps, ns = np.where(isP)[0], np.where(~isP)[0]
            pi, ni = ps[rng.integers(len(ps))], ns[rng.integers(len(ns))]
            old = E2*(1/D[pi, ps]).sum() + E2*(1/D[ni, ns]).sum() - E2/D[pi, pi] - E2/D[ni, ni]
            new = E2*(1/D[ni, ps]).sum() + E2*(1/D[pi, ns]).sum() - E2/D[ni, ni] - E2/D[pi, pi]
            if new < old or rng.random() < np.exp(-(new-old)/temp):
                isP[pi], isP[ni] = False, True
                Ec = coulomb(isP)
        bestEc = min(bestEc, Ec)
    return -(E_pair + T + bestEc)

# ============ 扩充网格: 偶偶 + 奇奇 + 奇A ============
GRID = []
# 奇奇核: A 偶, NZ 偶非零, Z=(A−NZ)/2 为奇 (A/2 与 NZ/2 奇偶不同)
for A in range(40, 241, 20):
    for NZ in (-14, -10, -6, -2, 2, 6, 10, 14):
        Z = (A-NZ)//2; N = (A+NZ)//2
        if Z >= 8 and N >= 8 and (Z % 2 == 1):  # 确保是奇奇
            GRID.append((Z, N))
# 奇A核: A 奇, NZ 奇
for A in range(41, 241, 20):
    for NZ in (-13, -9, -5, -1, 1, 5, 9, 13):
        Z = (A-NZ)//2; N = (A+NZ)//2
        if Z >= 8 and N >= 8 and A == Z+N:
            GRID.append((Z, N))
# 偶偶核: 复用原网格采样密度的一半 (A 步长 40) 控制总时长
for A in range(40, 241, 40):
    for NZ in (0, 4, 8, 16, 24):
        Z = (A-NZ)//2; N = (A+NZ)//2
        if Z >= 8 and N >= 8:
            GRID.append((Z, N))

GRID = sorted(set(GRID))
print(f"扩充网格: {len(GRID)} 点 (奇奇 {sum(1 for z,n in GRID if z%2==1 and (z+n)%2==0)}, "
      f"奇A {sum(1 for z,n in GRID if (z+n)%2==1)}, 偶偶 {sum(1 for z,n in GRID if z%2==0 and n%2==0)})")

Bv = []
for i, (Z, N) in enumerate(GRID):
    Bv.append(binding(Z, N))
    if (i+1) % 20 == 0:
        print(f"  ... {i+1}/{len(GRID)}")
Bv = np.array(Bv, float)
okm = ~np.isnan(Bv)
G = [g for g, m in zip(GRID, okm) if m]
B = Bv[okm]
Aarr = np.array([z+n for z, n in G], float)
Zarr = np.array([z for z, _ in G], float)
NZarr = np.array([n-z for z, n in G], float)

np.save('results/pairing_grid_B_v2.npy', B)
np.save('results/pairing_grid_AZN_v2.npy', np.array([[z+n, z, n] for z, n in G]))

X4 = np.vstack([Aarr, Aarr**(2/3), Zarr**2/Aarr**(1/3), NZarr**2/Aarr]).T
# 奇偶哑变量: 偶偶 +1, 奇A 0, 奇奇 -1
par = np.where(Aarr % 2 == 1, 0.0, np.where(Zarr % 2 == 0, 1.0, -1.0))
X5 = np.hstack([X4, par[:, None]])

c4, *_ = np.linalg.lstsq(X4, B, rcond=None)
c5, *_ = np.linalg.lstsq(X5, B, rcond=None)
r4, r5 = X4@c4 - B, X5@c5 - B

print()
print("="*84)
print("配对项提取结果 v2 (模型自有能量, 零实验借用)")
print("="*84)
print(f"四参拟合 RMS = {np.sqrt((r4**2).mean()):.2f} MeV")
print(f"五参(含奇偶) RMS = {np.sqrt((r5**2).mean()):.2f} MeV")
print(f"模型配对系数 delta_model = {c5[-1]:+.3f} MeV  (实验量级 ~12/sqrt(A)≈0.8–1.9)")
print()
print("四参拟合残差按奇偶分组:")
means = {}
for name, m in [('偶偶', (Aarr%2==0)&(Zarr%2==0)), ('奇A', Aarr%2==1), ('奇奇', (Aarr%2==0)&(Zarr%2==1))]:
    if m.any():
        means[name] = (float(r4[m].mean()), float(r4[m].std()), int(m.sum()))
        print(f"  {name}: n={m.sum():3d}  均值={r4[m].mean():+8.2f}  σ={r4[m].std():.2f}  均值标准误={r4[m].std()/np.sqrt(m.sum()):.2f}")
if len(means) == 3:
    mv = [means[nm][0] for nm in ('偶偶','奇A','奇奇')]
    spread = max(mv)-min(mv)
    se = max(means[nm][2] for nm in means)**-0.5 * np.sqrt((r4**2).mean())
    print(f"\n三组均值极差 = {spread:.2f} MeV; 最大的组均值标准误 ≈ {max(means[nm][1]/np.sqrt(means[nm][2]) for nm in means):.2f} MeV")
    verdict = '检测到奇偶交错信号' if spread > 3*max(means[nm][1]/np.sqrt(means[nm][2]) for nm in means) else '未检测到显著奇偶信号 (格点噪声淹没, 配对项无法从模型自有能量零借用提取)'
    print(f"裁决 → **{verdict}**")
json.dump(dict(n=int(len(B)), rms4=float(np.sqrt((r4**2).mean())), rms5=float(np.sqrt((r5**2).mean())),
               delta_model=float(c5[-1]), parity_means=means, verdict=verdict),
          open('results/pairing_probe_v2.json', 'w'), indent=1, ensure_ascii=False)
print("存档: results/pairing_probe_v2.json  pairing_grid_B_v2.npy  pairing_grid_AZN_v2.npy")
