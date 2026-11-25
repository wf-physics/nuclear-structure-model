# -*- coding: utf-8 -*-
"""
差分公式族 v4 — 模型自有液滴系数提取 (F0 公式闭合)
路线: 离散四键网络 B_model(A,Z) (v2 光滑格点引擎) → 网格采样 → 最小二乘提取
     B_model(A,Z) = a_v·A − a_s·A^(2/3) − a_c·Z²/A^(1/3) − a_sym·(N−Z)²/A
     四个系数全部由模型自身能量导出 (零实验借用) → 这就是模型的 F0 闭式公式
     格点壳噪声 (σ≈5 MeV, 离散化伪影非物理) 在拟合中平均掉
然后 F1-F4 差分公式全部用此闭式解析计算, 对照 AME 实验 → 适用边界
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
    pts = []
    for i in range(-m, m+1):
        for j in range(-m, m+1):
            for kk in range(-m, m+1):
                for f in fcc:
                    for b in basis:
                        pt = ((i+f[0]+b[0])*a, (j+f[1]+b[1])*a, (kk+f[2]+b[2])*a)
                        if np.linalg.norm(pt) <= R_c*1.02:
                            pts.append(pt)
    P = np.array(pts)
    if len(P) < A:
        for fac in (1.08, 1.15, 1.25, 1.40):
            extra = []
            for i in range(-m, m+1):
                for j in range(-m, m+1):
                    for kk in range(-m, m+1):
                        for f in fcc:
                            for b in basis:
                                pt = ((i+f[0]+b[0])*a, (j+f[1]+b[1])*a, (kk+f[2]+b[2])*a)
                                if np.linalg.norm(pt) <= R_c*fac:
                                    extra.append(pt)
            P = np.array(extra)
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

# ================= 第一步: 网格采样模型自有 B =================
print("="*84)
print("第一步: 离散四键网络 B_model 网格采样 (提取模型自有液滴系数)")
print("="*84)
GRID = []
for A in range(40, 241, 20):
    for dNZ in (0, 4, 8, 16, 24):
        for sgn in ((1,) if dNZ == 0 else (1, -1)):
            NZ = dNZ*sgn
            Z = (A - NZ)//2
            N = (A + NZ)//2
            if Z >= 8 and N >= 8 and A == Z+N:
                GRID.append((Z, N))
print(f"  采样点: {len(GRID)} 个 (A=40..240, N−Z ∈ [−24,24])")
Bv = []
for Z, N in GRID:
    b = binding(Z, N)
    Bv.append(b)
Bv = np.array(Bv, float)
okm = ~np.isnan(Bv)
GRID_ok = [g for g, m in zip(GRID, okm) if m]
Bv = Bv[okm]
print(f"  有效: {len(Bv)}")

# ================= 第二步: 模型自有液滴系数最小二乘 =================
Aarr = np.array([z+n for z, n in GRID_ok], float)
Zarr = np.array([z for z, _ in GRID_ok], float)
NZarr = np.array([n-z for z, n in GRID_ok], float)
X = np.vstack([Aarr, Aarr**(2/3), Zarr**2/Aarr**(1/3), NZarr**2/Aarr]).T
coef, res, *_ = np.linalg.lstsq(X, Bv, rcond=None)
a_v, a_s, a_c, a_sym = coef[0], -coef[1], -coef[2], -coef[3]
pred = X@coef
rms_fit = np.sqrt(np.mean((pred-Bv)**2))
print()
print("="*84)
print("第二步: 模型自有液滴系数 (全部由键网络能量导出, 零实验借用)")
print("="*84)
print(f"  B_model(A,Z) = {a_v:.3f}·A − {a_s:.3f}·A^(2/3) − {a_c:.4f}·Z²/A^(1/3) − {a_sym:.2f}·(N−Z)²/A")
print(f"  拟合 RMS = {rms_fit:.2f} MeV (即格点壳噪声被吸收后的残差)")
print(f"  对照: a_v 经验 15.8 | a_s 经验 18 (模型 {a_s:.1f}) | a_c 经验 0.71 (模型 {a_c:.3f}) | a_sym 经验 23 (模型 {a_sym:.1f})")
np.save('results/model_drop_coef.npy', coef)

def Bf(Z, N):
    A = Z+N
    return (a_v*A + coef[1]*A**(2/3) + coef[2]*Z**2/A**(1/3) + coef[3]*(N-Z)**2/A)

def Bexp(Z, N):
    v = AME.get((Z, N))
    return v['BEA']*v['A']/1000.0 if v and v.get('BEA') else None

def judge(dev, Qe):
    return ('过' if abs(dev) <= max(0.4, 0.2*abs(Qe)) else
            '半' if abs(dev) <= max(1.0, 0.5*abs(Qe)) else '不过')

rows = {'F1_Sn': [], 'F2_Qbm': [], 'F3_Qa': [], 'F4_Qfs': []}

# ================= F1-F4 闭式差分对照 =================
print()
print("="*84)
print("第三步: F1-F4 差分公式 (闭式, 零噪声) 对照 AME 实验 → 适用边界")
print("="*84)

print("--- F1  S_n = B(A,Z) − B(A−1,Z) ---")
for Z, N, nm in [(20,20,'Ca-40'),(20,28,'Ca-48'),(26,30,'Fe-56'),(28,30,'Ni-58'),
                 (30,34,'Zn-64'),(36,44,'Kr-80'),(38,50,'Sr-88'),(40,50,'Zr-90'),
                 (50,62,'Sn-112'),(82,126,'Pb-208')]:
    v1, v2 = AME.get((Z,N)), AME.get((Z,N-1))
    if not (v1 and v2 and v1.get('BEA') and v2.get('BEA')):
        continue
    Sn_m = Bf(Z, N) - Bf(Z, N-1)
    Sn_e = v1['BEA']*v1['A']/1000. - v2['BEA']*v2['A']/1000.
    dev = Sn_m - Sn_e
    vd = '过' if abs(dev) <= 2.0 else ('半' if abs(dev) <= 5.0 else '不过')
    print(f"  {nm:>8s} 模型={Sn_m:8.2f} 实验={Sn_e:8.2f} 偏差={dev:+8.2f} {vd}")
    rows['F1_Sn'].append(dict(nucl=nm, Sn_m=Sn_m, Sn_e=Sn_e, dev=dev, v=vd))

print("--- F2  Q_β⁻ = B(A,Z+1) − B(A,Z) − 0.782 ---")
for Z, N, nm in [(6,8,'C-14'),(11,13,'Na-24'),(15,17,'P-32'),(16,19,'S-35'),
                 (19,23,'K-42'),(25,31,'Mn-56'),(26,33,'Fe-59'),
                 (38,52,'Sr-90'),(55,82,'Cs-137'),(61,86,'Pm-147')]:
    e = [AME.get(kk) for kk in [(Z,N),(Z+1,N-1)]]
    if any(v is None or v.get('ME') is None for v in e):
        continue
    Qm = Bf(Z+1, N-1) - Bf(Z, N) - 0.782   # β⁻: (A,Z)→(A,Z+1), A 不变
    Qe = (e[0]['ME']-e[1]['ME'])/1000.
    dev = Qm - Qe
    vd = judge(dev, Qe)
    print(f"  {nm:>8s} 模型={Qm:8.3f} 实验={Qe:8.3f} 偏差={dev:+8.3f} {vd}")
    rows['F2_Qbm'].append(dict(nucl=nm, Qm=Qm, Qe=Qe, dev=dev, v=vd))

print("--- F3  Q_α = B(A−4,Z−2) + B_α(锚28.296) − B(A,Z) ---")
for Z, N, nm in [(62,85,'Sm-147'),(64,88,'Gd-152'),(74,106,'W-180'),(84,126,'Po-210'),
                 (86,136,'Rn-222'),(88,138,'Ra-226'),(90,142,'Th-232'),
                 (92,143,'U-235'),(92,146,'U-238'),(94,145,'Pu-239'),
                 (94,146,'Pu-240'),(95,146,'Am-241'),(96,148,'Cm-244'),(98,154,'Cf-252')]:
    v1, v2 = AME.get((Z,N)), AME.get((Z-2,N-2))
    if not (v1 and v2 and v1.get('ME') is not None and v2.get('ME') is not None):
        continue
    Qm = Bf(Z-2, N-2) + 28.2957 - Bf(Z, N)
    Qe = (v1['ME']-v2['ME']-2424.916)/1000.
    dev = Qm - Qe
    vd = judge(dev, Qe)
    print(f"  {nm:>8s} 模型={Qm:8.2f} 实验={Qe:8.2f} 偏差={dev:+8.2f} {vd}")
    rows['F3_Qa'].append(dict(nucl=nm, Qm=Qm, Qe=Qe, dev=dev, v=vd))

print("--- F4  Q_fs = B(f1)+B(f2) − B(母核) [BEA 口径] ---")
FISSIONS = [
    (92,146,'U-238', [(54,86,'Xe-140'),(38,57,'Sr-95')], 3),
    (92,146,'U-238', [(46,73,'Pd-119'),(46,73,'Pd-119')], 0),
    (92,143,'U-235', [(56,85,'Ba-141'),(36,53,'Kr-89')], 3),
    (94,146,'Pu-240', [(57,88,'La-145'),(37,55,'Rb-92')], 3),
    (98,154,'Cf-252', [(56,86,'Ba-142'),(42,65,'Mo-107')], 3),
]
for Z, N, nm, frags, xn in FISSIONS:
    vp = AME.get((Z, N))
    if not (vp and vp.get('BEA')):
        continue
    Qm = -Bf(Z, N) + sum(Bf(fz, fn) for fz, fn, _ in frags)
    Qe = -vp['BEA']*vp['A']/1000.
    ok = True
    for fz, fn, _ in frags:
        vf = AME.get((fz, fn))
        if vf is None or vf.get('BEA') is None:
            ok = False; break
        Qe += vf['BEA']*vf['A']/1000.
    if not ok:
        print(f"  {nm:>8s} 碎片实验质量缺"); continue
    dev = Qm - Qe
    vd = judge(dev, Qe)
    frag_s = '+'.join(f[2] for f in frags) + (f'+{xn}n' if xn else '')
    print(f"  {nm:>8s} {frag_s:>18s} 模型={Qm:8.1f} 实验={Qe:8.1f} 偏差={dev:+8.1f} {vd}")
    rows['F4_Qfs'].append(dict(parent=nm, frags=frag_s, Qm=Qm, Qe=Qe, dev=dev, v=vd))

print()
print("="*84)
print("汇总: 差分公式族适用边界 (v4, 模型自有闭式)")
print("="*84)
for fid in ['F1_Sn', 'F2_Qbm', 'F3_Qa', 'F4_Qfs']:
    t = rows[fid]
    if not t:
        continue
    ok = sum(1 for r in t if r['v'] == '过')
    hf = sum(1 for r in t if r['v'] == '半')
    mae = np.mean([abs(r['dev']) for r in t])
    print(f"  {fid:8s}: n={len(t):2d}  过={ok}  半={hf}  不过={len(t)-ok-hf}   MAE={mae:.2f} MeV")

json.dump(rows, open('results/q_formula_v4.json', 'w'), indent=1, ensure_ascii=False)
json.dump(dict(a_v=float(a_v), a_s=float(a_s), a_c=float(a_c), a_sym=float(a_sym),
               fit_rms=float(rms_fit), n_grid=len(Bv)),
          open('results/model_drop_coef.json', 'w'), indent=1)
print("存档: results/q_formula_v4.json  model_drop_coef.{npy,json}")
