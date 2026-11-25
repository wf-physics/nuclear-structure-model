# -*- coding: utf-8 -*-
"""
键能网络轻核标定扫描（v4.0 §6 第一阶段的量化预演）
模型（v4.0 §1.1 原样）:
    E_bond(r) = 5·(2.2 - r) MeV（即 7·(2.2-r)/1.4）
    S(k) = min(1, 2/(k+0.1))
    E_eff = S(k)·P_ij·E_bond(r),  P = 1 (pn) 或 w_pp (pp/nn)
    结合能 = B_total − E_Coul（³He、⁴He 含 pp 库仑排斥 1.44/r）

标定对象（AME2020 结合能）:
    ²H 2.2245 / ³H 8.4818 / ³He 7.7180 / ⁴He 28.2957 MeV
几何假设：均匀键长 r（²H 一条边；³H/³He 三角形 3 条 pn 边，k=2；
          ⁴He 四面体 6 条边 = 4 pn + 2 pp + 2 nn，k=3）
"""
import numpy as np
import itertools, json

EXP = {"2H": 2.2245, "3H": 8.4818, "3He": 7.7180, "4He": 28.2957}
S2, S3 = 2/2.1, 2/3.1   # k=2, k=3 的软饱和因子

def model_B(nuc, r, w_pp):
    eb = 5.0 * (2.2 - r)
    if nuc == "2H":  return eb * 1.0                       # k=1, S=1
    if nuc == "3H":  return 3 * eb * S2                    # 3 pn, k=2
    if nuc == "3He": return 3 * eb * S2 - 1.44 / r         # + pp 库仑排斥
    if nuc == "4He":
        B = 4 * eb * S3 + 4 * eb * S3 * w_pp               # 4 pn + 2pp+2nn, k=3
        return B - 1.44 / r                                # + pp 库仑排斥
    raise ValueError(nuc)

def nrmse(residuals, refs):
    res = np.array(list(residuals.values())) if isinstance(residuals, dict) else np.array(residuals)
    ref = np.array(refs)
    return float(np.sqrt(np.mean(res**2)) / np.sqrt(np.mean(ref**2)))

def main():
    r_grid  = np.arange(0.80, 2.2001, 0.005)
    w_grid  = np.arange(0.60, 0.8501, 0.005)
    best = None
    for r, w in itertools.product(r_grid, w_grid):
        preds = {k: model_B(k, r, w) for k in EXP}
        res = {k: preds[k] - EXP[k] for k in EXP}
        score = nrmse(res, list(EXP.values()))
        if best is None or score < best[0]:
            best = (score, r, w, preds, res)

    score, r, w, preds, res = best
    print("=" * 86)
    print("键能网络轻核标定：全局最优 (r, w_pp) 扫描结果")
    print("=" * 86)
    print(f"最优: r = {r:.3f} fm, w_pp = {w:.3f}, NRMSE = {score*100:.1f}%\n")
    print(f"{'核素':<8}{'实验 (MeV)':>12}{'模型 (MeV)':>12}{'残差 (MeV)':>12}")
    for k in EXP:
        print(f"{k:<8}{EXP[k]:>12.4f}{preds[k]:>12.4f}{res[k]:>12.4f}")

    print("\n--- 各核素'单独最优键长'（w_pp 扫描取最优）---")
    print(f"{'核素':<8}{'最优 r (fm)':>12}{'最优 w_pp':>10}{'残差 (MeV)':>12}")
    for k in EXP:
        bb = None
        for r, w in itertools.product(r_grid, w_grid):
            d = model_B(k, r, w) - EXP[k]
            if bb is None or abs(d) < abs(bb[3]):
                bb = (r, w, model_B(k, r, w), d)
        print(f"{k:<8}{bb[0]:>12.3f}{bb[1]:>10.3f}{bb[3]:>12.4f}")

    print("""
结论摘要（对 v4.0 §6 第一阶段的意义）：
1. 单一 (r, w_pp) 无法同时拟合 ²H/³H/³He/⁴He：²H 想要 r≈1.8 fm，
   ³H/³He 想要 r≈1.6 fm，而 ⁴He 需要 r≈1.0 fm 才凑够 28.3 MeV ——
   软饱和因子 S(3)=0.645 把 ⁴He 每条边压到 ³H 的 ~68%，与实验趋势相反
   （实验上 ⁴He 每核子结合能 7.07 MeV ≫ ³H 的 2.83 MeV）。
2. 这不是参数没调好，而是结构性冲突：纯距离两体可加 + 饱和折扣，
   天然给出'越密排总键能越低'的走向；真实核结合能来自饱和体积项 +
   多体/动量相关贡献。该冲突触发了 §8 证伪条件 1 的边缘情形。
3. 可能的出路：(i) 允许 r 随 A 变化（放弃统一键长）；
   (ii) S(k) 反向（配位多→折扣小）；(iii) 加体积项 B ∝ A。
   建议在 v4.1 讨论段明确写出其中一条并重新走 §6 三阶段标定。
""")

if __name__ == "__main__":
    main()
