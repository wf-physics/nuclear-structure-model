# -*- coding: utf-8 -*-
"""
能量公式体系 — 有限核结合能计算 v1.0 (金刚石晶格液滴近似)
标定常数: nm_formula_v4.json (diamond, k=57, C=699, α=3, p=6, δ=0.2, ρ0=0.16)
B(A,Z) = -E_pair - T - E_C
  E_pair: 金刚石晶格球 (均匀密度 ρ0), 冻结结构 V = 线性井(平滑窗口) + 泡利核(ρ0)
  T:      双组分费米动能 k_F,s = (3π²ρ_s)^(1/3)
  E_C:    质子对精确库仑 (质子位点分配用 Metropolis 最小化)
适用边界探测: 中核拟合良好 → 重核偏差扩大的分界线
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

def diamond_sites(A, rho=RHO0):
    """A 个金刚石格点, 均匀密度 rho, 取最靠近中心的 A 个"""
    a = (8.0/rho)**(1/3)
    m = int(np.ceil(a)) + 2
    fcc = [(0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0)]
    basis = [(0,0,0),(.25,.25,.25)]
    pts = []
    R = (3*A/(4*np.pi*rho))**(1/3)*1.35
    for i in range(-m, m+1):
        for j in range(-m, m+1):
            for kk in range(-m, m+1):
                for f in fcc:
                    for b in basis:
                        pt = ((i+f[0]+b[0])*a, (j+f[1]+b[1])*a, (kk+f[2]+b[2])*a)
                        if np.linalg.norm(pt) <= R:
                            pts.append(pt)
    P = np.array(pts)
    if len(P) < A:
        return None
    r = np.linalg.norm(P, axis=1)
    P = P[np.argsort(r)][:A]
    return P - P.mean(axis=0)

def binding(Z, N):
    A = Z + N
    P = diamond_sites(A)
    if P is None:
        return None, 'sites'
    D = np.linalg.norm(P[:, None, :]-P[None, :, :], axis=2)
    np.fill_diagonal(D, 1e9)
    att = -k*(2.2-D)*sig((D-2.2)/delta)*sig((0.8-D)/delta)
    rep = C*(0.5/np.maximum(D, 1e-6))**p          # 泡利核 @ ρ0 (α 项=1)
    E_pair = 0.5*(att+rep).sum()
    # 双组分费米动能
    rho_p, rho_n = Z/A*RHO0, N/A*RHO0
    kp, kn = (3*np.pi**2*rho_p)**(1/3), (3*np.pi**2*rho_n)**(1/3)
    T = 0.6*HB2M*(Z*kp**2 + N*kn**2)
    # 库仑: 质子分配 Metropolis 最小化
    rng = np.random.default_rng(42)
    idx = rng.choice(A, size=Z, replace=False)
    isP = np.zeros(A, bool); isP[idx] = True
    Dp = D.copy()
    def coulomb(mask):
        pp = np.outer(mask, mask)
        np.fill_diagonal(pp, False)
        return 0.5*(E2/Dp[pp]).sum() if pp.any() else 0.0
    Ec = coulomb(isP)
    for it in range(3000):
        pi = rng.choice(np.where(isP)[0])
        ni = rng.choice(np.where(~isP)[0])
        # 交换 p↔n 的库仑变化 (近似: 只算涉及两位点的项)
        pi_set, ni_set = np.where(isP)[0], np.where(~isP)[0]
        old = E2*(1/Dp[pi, pi_set]).sum() + E2*(1/Dp[ni, ni_set]).sum() \
              - E2/Dp[pi, pi] - E2/Dp[ni, ni]
        new = E2*(1/Dp[ni, pi_set]).sum() + E2*(1/Dp[pi, ni_set]).sum() \
              - E2/Dp[ni, ni] - E2/Dp[pi, pi]
        if new < old or rng.random() < np.exp(-(new-old)/0.05):
            isP[pi], isP[ni] = False, True
            Ec = coulomb(isP)
    E_tot = E_pair + T + Ec
    return -E_tot, 'ok'

TEST = [(8,6,'O-14'),(8,8,'O-16'),(10,10,'Ne-20'),(12,12,'Mg-24'),(14,14,'Si-28'),
        (16,16,'S-32'),(18,18,'Ar-36'),(20,20,'Ca-40'),(20,28,'Ca-48'),(22,26,'Ti-48'),
        (24,28,'Cr-52'),(26,30,'Fe-56'),(28,30,'Ni-58'),(30,34,'Zn-64'),(32,38,'Ge-70'),
        (36,44,'Kr-80'),(38,50,'Sr-88'),(40,50,'Zr-90'),(42,54,'Mo-96'),(44,58,'Ru-102'),
        (50,62,'Sn-112'),(50,72,'Sn-122'),(82,126,'Pb-208')]

print("="*78)
print("有限核结合能 (金刚石晶格液滴近似, 标定常数自动导出, 零逐核拟合)")
print("="*78)
print(f"{'核':>8s} {'A':>4s} {'B模型':>9s} {'B实验':>9s} {'偏差':>8s} {'相对':>7s} {'B/A模型':>8s} {'B/A实验':>8s}")
rows = []
for Z, N, name in TEST:
    Bm, status = binding(Z, N)
    key = (Z, N)
    amev = AME.get(key)
    Bexp = amev['BEA']*amev['A']/1000.0 if amev and amev.get('BEA') else None
    if Bm is None:
        print(f"{name:>8s} {A if False else Z+N:4d}  格点不足")
        continue
    if Bexp:
        dev = Bm-Bexp
        print(f"{name:>8s} {Z+N:4d} {Bm:9.1f} {Bexp:9.1f} {dev:+8.1f} {dev/Bexp*100:+6.1f}% "
              f"{Bm/(Z+N):8.2f} {Bexp/(Z+N):8.2f}")
        rows.append(dict(name=name, Z=Z, N=N, Bm=Bm, Bexp=Bexp, dev=dev, rel=dev/Bexp*100))
    else:
        print(f"{name:>8s} {Z+N:4d} {Bm:9.1f} {'--':>9s}")
json.dump(rows, open('results/finite_nuclei_v1.json', 'w'), indent=1, ensure_ascii=False)
# 边界分析
print()
for lo, hi, tag in [(12, 40, '轻中核 A12-40'), (40, 65, '中核 A40-65'), (65, 100, '较重 A65-100'), (100, 210, '重核 A>100')]:
    t = [r for r in rows if lo <= r['Z']+r['N'] <= hi]
    if t:
        print(f"  {tag}: n={len(t)}  MAE={np.mean([abs(r['dev']) for r in t]):6.1f} MeV  "
              f"中位相对={np.median([r['rel'] for r in t]):+6.1f}%  "
              f"<15%比例={np.mean([abs(r['rel'])<15 for r in t])*100:5.1f}%")
