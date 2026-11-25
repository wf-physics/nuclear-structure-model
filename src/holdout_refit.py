#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
六参数宏观-微观全局拟合 + 留出验证（Z>=88 全剔除，规则冻结）。

复现目标（预印本 §3.3）:
  B(Z,N) = aV*A - aS*A^(2/3) - aC*Z(Z-1)/A^(1/3) - aA*(N-Z)^2/A + aP*d/sqrt(A) + dB_shell
  dB_shell = C_shell*[ c*A^(1/3) - (F(N)+F(Z))/(A/2)^(2/3) ]
  F(X) = (3/5)[ (M_i^(5/3)-M_{i-1}^(5/3))/(M_i-M_{i-1})*(X-M_{i-1}) - (X^(5/3)-M_{i-1}^(5/3)) ]
  magic = {2,8,14,28,50,82,126,184}
  样本内: MAE 2.297 / RMS 3.069, aV=15.287 aS=16.307 aC=0.690 aA=22.251 aP=11.049 C=3.226

留出规则（冻结, 拟合前先写下）: 剔除 Z>=88 的全部核素, 在剩余核素上重拟 6 参数,
再计算 238U 的 B / S_n / Q_alpha / Q_fs(Pd+Pd) 四个偏差。
"""
import numpy as np, json, re, sys
from scipy.optimize import least_squares

DATA = "/workspace/核结构研究/repo/data/mass_1.mas20"
MAGIC = [2, 8, 14, 28, 50, 82, 126, 184]
B_ALPHA = 28.296  # MeV, alpha 粒子结合能锚点（论文唯一显式实验借用）


def parse_masses(path, only_experimental=True):
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    for ln in lines[15:]:  # 跳过文件头（MASS LIST 表头从 15 行附近开始）
        if len(ln) < 67:
            continue
        try:
            n = int(ln[4:9]); z = int(ln[9:14]); a = int(ln[14:19])
        except ValueError:
            continue
        if a == 0 or z == 0:
            continue
        bea_field = ln[54:67]
        bea_str = bea_field.strip().replace("#", "")
        try:
            bea = float(bea_str)
        except ValueError:
            continue
        if only_experimental and "#" in bea_field:
            continue
        be_total = bea * a / 1000.0  # MeV
        rows.append((z, n, a, be_total))
    return rows


def F_piecewise(X):
    """Myers-Swiatecki 平滑函数 F(X)，向量化。"""
    X = np.maximum(np.asarray(X, float), 1e-9)
    out = np.zeros_like(X)
    # 低于首个幻数：用 [0, M0] 区间外推（M_{i-1}=0）
    m0 = MAGIC[0]
    slope0 = (m0 ** (5 / 3) - 0.0) / m0
    low = X <= m0
    out[low] = (3 / 5) * (slope0 * X[low] - X[low] ** (5 / 3))
    for i in range(1, len(MAGIC)):
        Ml, Mr = MAGIC[i - 1], MAGIC[i]
        slope = (Mr ** (5 / 3) - Ml ** (5 / 3)) / (Mr - Ml)
        seg = (X > Ml) & (X <= Mr)
        out[seg] = (3 / 5) * (slope * (X[seg] - Ml) - (X[seg] ** (5 / 3) - Ml ** (5 / 3)))
    # 高于最后一个幻数：最后一个区间斜率外推
    Ml, Mr = MAGIC[-2], MAGIC[-1]
    slope = (Mr ** (5 / 3) - Ml ** (5 / 3)) / (Mr - Ml)
    high = X > Mr
    out[high] = (3 / 5) * (slope * (X[high] - Mr) - (X[high] ** (5 / 3) - Mr ** (5 / 3)))
    return out


def shell_corr(A, N, Z, C, c):
    return C * (c * A ** (1 / 3) - (F_piecewise(N) + F_piecewise(Z)) / (A / 2) ** (2 / 3))


def B_model(Z, N, p, c):
    aV, aS, aC, aA, aP, C = p
    A = Z + N
    ee = (Z % 2 == 0) & (N % 2 == 0)
    oo = (Z % 2 == 1) & (N % 2 == 1)
    d = np.where(ee, 1.0, np.where(oo, -1.0, 0.0))
    return (aV * A - aS * A ** (2 / 3) - aC * Z * (Z - 1) / A ** (1 / 3)
            - aA * (N - Z) ** 2 / A + aP * d / np.sqrt(A)
            + shell_corr(A, N, Z, C, c))


def fit(rows, c, p0=None):
    Z = np.array([r[0] for r in rows], float)
    N = np.array([r[1] for r in rows], float)
    Bexp = np.array([r[3] for r in rows])
    if p0 is None:
        p0 = np.array([15.5, 17.0, 0.7, 23.0, 11.0, 3.0])
    def resid(p):
        return B_model(Z, N, p, c) - Bexp
    sol = least_squares(resid, p0, method="lm", max_nfev=20000)
    return sol.x, resid(sol.x)


def stats(resid):
    return float(np.mean(np.abs(resid))), float(np.sqrt(np.mean(resid ** 2)))


def main():
    rows_all = parse_masses(DATA, only_experimental=False)
    rows16 = [r for r in rows_all if r[2] >= 16]
    n_exp = len([r for r in rows16 if True])
    print(f"nuclei total={len(rows_all)}, A>=16: {len(rows16)} (paper says 3487)")

    results = {}
    for c in [1.0, 0.5, 0.866, 0.417, 2.0]:
        p, r = fit(rows16, c)
        mae, rms = stats(r)
        print(f"c={c}: params={np.round(p,4)}  MAE={mae:.3f} RMS={rms:.3f}")
        results[c] = p

    # 也尝试包含估值核素
    rows_all_est = parse_masses(DATA, only_experimental=False)
    rows16_est = [r for r in rows_all_est if r[2] >= 16]
    print(f"\nwith estimated included: A>=16 count={len(rows16_est)}")
    for c in [1.0, 0.5]:
        p, r = fit(rows16_est, c)
        mae, rms = stats(r)
        print(f"c={c}: params={np.round(p,4)}  MAE={mae:.3f} RMS={rms:.3f}")


if __name__ == "__main__":
    main()
