# -*- coding: utf-8 -*-
"""
β 衰变相空间因子 f 的精确数值计算：自由中子 vs 氚(³H)
目的：定量检验《核结构拓扑网络理论法则 v4.0》"事例2"中的相空间论证
谱形（允许型跃迁）: dN/dT ∝ F(Z_f, E) · p · E · (Q - T)^2
  T  = 电子动能 (MeV)
  E  = T + m_e（总能量）
  p  = sqrt(T(T+2m_e))（动量，MeV/c）
  F  = 费米函数（非相对论近似）：F = 2πη / (1 - exp(-2πη))，η = Z_f·α·E/p
       β⁻（电子受子核吸引）时 F > 1 —— 这正是"库仑吸引增强电子波函数在核处权重"的标准表述
"""
import numpy as np

M_E = 0.510998950        # MeV, 电子质量
ALPHA = 1.0 / 137.035999 # 精细结构常数

def fermi_nonrel(Z, T):
    """非相对论费米函数 F(Z, T)"""
    T = np.asarray(T, dtype=float)
    p = np.sqrt(T * (T + 2 * M_E))
    W = T + M_E
    eta = Z * ALPHA * W / np.maximum(p, 1e-12)
    x = 2 * np.pi * eta
    # x 大时 F ≈ x（极限解析形式），避免数值问题
    F = np.where(x > 700, x, x / (-np.expm1(-x)))
    return F

def f_value(Q, Z, use_fermi=True, n=400001):
    """相空间积分 f = ∫ F·p·W·(Q-T)^2 dT （同一口径下比较，量纲不影响比值）"""
    T = np.linspace(0.0, Q, n)
    T[0] = 1e-12
    p = np.sqrt(T * (T + 2 * M_E))
    W = T + M_E
    base = p * W * (Q - T) ** 2
    if use_fermi:
        y = base * fermi_nonrel(Z, T)
    else:
        y = base
    try:
        val = float(np.trapezoid(y, T))
    except AttributeError:
        val = float(np.trapz(y, T))
    # 换算到标准无量纲约定 f_std = ∫ F p_nat W_nat (W0-W)_nat² dW_nat
    # （每 MeV 量除以 m_e，共 5 个因子：p, W, Δ², dT）→ 可与文献 log ft 直接对比
    return val / M_E**5

def main():
    Q_n = 0.782343   # MeV, 自由中子 β⁻ Q 值
    Q_t = 0.018591   # MeV, 氚 β⁻ Q 值
    Z_n = 1          # 子核电荷：质子
    Z_t = 2          # 子核电荷：³He

    tau_n = 879.4              # s, 自由中子平均寿命（实验）
    t12_n = tau_n * np.log(2)  # s, 半衰期 ≈ 609.6 s
    t12_t = 12.32 * 3.1557e7   # s, 氚半衰期 = 3.888e8 s

    f_n = f_value(Q_n, Z_n, use_fermi=True)
    f_t = f_value(Q_t, Z_t, use_fermi=True)
    f_n0 = f_value(Q_n, Z_n, use_fermi=False)
    f_t0 = f_value(Q_t, Z_t, use_fermi=False)

    lam_ratio_exp = (1.0 / tau_n) / (np.log(2) / t12_t)  # λ_n / λ_t 实验值
    q5_naive = (Q_n / Q_t) ** 5

    print("=" * 78)
    print("β 相空间因子精确积分（检验 v4.0 事例2 的相空间论证）")
    print("=" * 78)
    print(f"{'量':<46}{'数值':>16}")
    print("-" * 78)
    print(f"f(自由中子), 含费米函数 F(Z=1)      {f_n:>14.4f}")
    print(f"f(氚 ³H),   含费米函数 F(Z=2)      {f_t:>14.4e}")
    print(f"f(自由中子), 裸谱 F=1               {f_n0:>14.4f}")
    print(f"f(氚 ³H),   裸谱 F=1               {f_t0:>14.4e}")
    print("-" * 78)
    print(f"粗略标度 (Q_n/Q_t)^5（纲领中的 1.3e8）  {q5_naive:>14.4e}")
    print(f"精确相空间比 f_n/f_t（裸谱）            {f_n0/f_t0:>14.4e}")
    print(f"精确相空间比 f_n/f_t（含费米函数）      {f_n/f_t:>14.4e}")
    print(f"实验寿命比 λ_n/λ_t（用平均寿命口径）    {lam_ratio_exp:>14.4e}")
    print("-" * 78)
    ft_n = f_n * t12_n
    ft_t = f_t * t12_t
    print(f"log(ft) 自由中子  = {np.log10(ft_n):.3f}   （文献值 ≈ 3.01）")
    print(f"log(ft) 氚 ³H     = {np.log10(ft_t):.3f}   （文献值 ≈ 3.05）")
    print(f"若矩阵元相等，预言 λ 比 = {f_n/f_t:.3e}；实验 = {lam_ratio_exp:.3e}")
    resid = (f_n/f_t) / lam_ratio_exp
    print(f"残差因子 = {resid:.3f}  →  相空间+库仑(费米函数)已解释寿命比的"
          f"绝大部分，残差 {abs(resid-1)*100:.0f}% 属矩阵元/耦合常数细节")

if __name__ == "__main__":
    main()
