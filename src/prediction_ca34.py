#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""预言卡 ⁹C+²⁵Si→³⁴Ca+γ：一条龙复现脚本。

流程:
  1. 从 AME2020 原文件解析实验/估值核素 (列位: N=4:9, Z=9:14, A=14:19, B/A=54:67);
  2. 仅实验核素 (A>=16, 2484 核) 拟合六参数宏观-微观模型 (c=0.25 固定);
  3. 局部误差棒: Z=16-24, N=8-20 实测核素残差 RMS;
  4. 输出预言数字: 聚变 Q / 2p 衰变能 / S_p / B(34Ca);
  5. 对照 NUBASE2020 束缚核过滤器 (repo/results/nubase_bound.json) 的通道枚举
     见 fusion_channels_bound.json 生成段 (本脚本输出主数字)。

输出: results/prediction_ca34_numbers.json
论文引用: 预言卡_9C25Si_to_34Ca_2026-11-24.md
"""
import json, os
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data", "mass_1.mas20")
MAGIC = [2, 8, 14, 28, 50, 82, 126, 184]
B_ALPHA = 28.296
C_FIXED = 0.25


def parse(path):
    exp, est = {}, {}
    for ln in open(path, encoding="utf-8", errors="replace").readlines()[15:]:
        if len(ln) < 67:
            continue
        try:
            n = int(ln[4:9]); z = int(ln[9:14]); a = int(ln[14:19])
        except ValueError:
            continue
        if a == 0 or z == 0:
            continue
        f = ln[54:67]
        try:
            bea = float(f.strip().replace("#", ""))
        except ValueError:
            continue
        (est if "#" in f else exp)[(z, n)] = bea * (z + n) / 1000.0
    return exp, est


def F_piecewise(X):
    X = np.maximum(np.asarray(X, float), 1e-9)
    out = np.zeros_like(X)
    m0 = MAGIC[0]
    slope0 = m0 ** (5 / 3) / m0
    low = X <= m0
    out[low] = (3 / 5) * (slope0 * X[low] - X[low] ** (5 / 3))
    for i in range(1, len(MAGIC)):
        Ml, Mr = MAGIC[i - 1], MAGIC[i]
        slope = (Mr ** (5 / 3) - Ml ** (5 / 3)) / (Mr - Ml)
        seg = (X > Ml) & (X <= Mr)
        out[seg] = (3 / 5) * (slope * (X[seg] - Ml) - (X[seg] ** (5 / 3) - Ml ** (5 / 3)))
    Ml, Mr = MAGIC[-2], MAGIC[-1]
    slope = (Mr ** (5 / 3) - Ml ** (5 / 3)) / (Mr - Ml)
    high = X > Mr
    out[high] = (3 / 5) * (slope * (X[high] - Mr) - (X[high] ** (5 / 3) - Mr ** (5 / 3)))
    return out


def B_model(Z, N, p):
    aV, aS, aC, aA, aP, C = p
    A = Z + N
    Z = np.asarray(Z, float); N = np.asarray(N, float); A = np.asarray(A, float)
    ee = (Z % 2 == 0) & (N % 2 == 0)
    oo = (Z % 2 == 1) & (N % 2 == 1)
    d = np.where(ee, 1.0, np.where(oo, -1.0, 0.0))
    return (aV * A - aS * A ** (2 / 3) - aC * Z * (Z - 1) / A ** (1 / 3)
            - aA * (N - Z) ** 2 / A + aP * d / np.sqrt(A)
            + C * (C_FIXED * A ** (1 / 3) - (F_piecewise(N) + F_piecewise(Z)) / (A / 2) ** (2 / 3)))


def main():
    exp, est = parse(DATA)
    rows = [(z, n, z + n, b) for (z, n), b in exp.items() if z + n >= 16]
    Z = np.array([r[0] for r in rows], float)
    N = np.array([r[1] for r in rows], float)
    Bexp = np.array([r[3] for r in rows])
    sol = least_squares(lambda p: B_model(Z, N, p) - Bexp,
                        np.array([15.4, 16.8, 0.69, 22.4, 11.7, 3.5]),
                        method="lm", max_nfev=20000)
    p = sol.x
    r = B_model(Z, N, p) - Bexp
    sig_local = float(np.sqrt(np.mean(r[(Z >= 16) & (Z <= 24) & (N >= 8) & (N <= 20)] ** 2)))

    def bp(z, n):
        return float(B_model(np.array([z], float), np.array([n], float), p)[0])

    B34 = bp(20, 14)
    out = {
        "model_params_expfit": {k: round(float(v), 4) for k, v in
                                zip(["aV", "aS", "aC", "aA", "aP", "C_shell"], p)},
        "n_fit_experimental": len(rows),
        "sigma_local_Z16_24_N8_20": round(sig_local, 3),
        "B34Ca_pred_MeV": round(B34, 2),
        "fusion_Q_9C_25Si": round(B34 - exp[(6, 3)] - exp[(14, 11)], 2),
        "fusion_Q_10C_24Si": round(B34 - exp[(6, 4)] - exp[(14, 10)], 2),
        "S2p_ET_34Ca_to_32Ar": round(exp[(18, 14)] - B34, 2),
        "Sp_model_34Ca_to_33K": round(bp(19, 14) - B34, 2),
    }
    json.dump(out, open(os.path.join(HERE, "..", "results", "prediction_ca34_numbers.json"), "w"),
              indent=1, ensure_ascii=False)
    print(json.dumps(out, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
