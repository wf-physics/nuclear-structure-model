# -*- coding: utf-8 -*-
"""引擎核心 (复制自 q_formula_v4.py, 独立导入用, 未改动任何常数)"""
import numpy as np, json
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

